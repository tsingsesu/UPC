/* Standalone isotropic Gaussian KDE; no solver, DOM or external dependency.
 * Keep the complete factory serializable via createDensity.toString() for Workers.
 * A KDE describes declared simulation samples; it is not experimentally verified
 * football probability. The bandwidth is a smoothing scale in metres.
 */
'use strict';
function createDensity() {
  const VERSION = 'upc-density-1.1.0';
  const SQRT_2PI = Math.sqrt(2 * Math.PI);
  const finite = Number.isFinite;
  const failure = (message) => ({version:VERSION,status:'invalid_input',issues:[message]});
  function plainObject(value) { return value !== null && typeof value === 'object' && !Array.isArray(value); }

  // For |x|<3, integrate the Gaussian with a convergent positive-term series.
  // For x>=3, use the Laplace continued fraction for Q(x)/phi(x).
  // No fitted third-party polynomial coefficients are embedded. Returning tails
  // directly avoids 1-CDF cancellation in small, far-tail rectangles.
  function positiveNormalTail(x) {
    if (x === Infinity || x > 38.5) return 0;
    const phi = Math.exp(-0.5*x*x)/SQRT_2PI;
    if (x >= 3) {
      let fraction = 0;
      for (let k = 128; k >= 1; k--) fraction = k/(x+fraction);
      return phi/(x+fraction);
    }
    let term = x, sum = x;
    for (let k = 1; k < 256; k++) {
      term *= x*x/(2*k+1);
      sum += term;
      if (term <= Number.EPSILON * Math.abs(sum)) break;
    }
    return Math.max(0, Math.min(0.5, 0.5-phi*sum));
  }
  function normalCDF(x) {
    if (typeof x !== 'number' || Number.isNaN(x)) return NaN;
    return x >= 0 ? 1-positiveNormalTail(x) : positiveNormalTail(-x);
  }
  function normalInterval(low, high) {
    if (low >= high) return 0;
    let result;
    if (low >= 0) result = positiveNormalTail(low)-positiveNormalTail(high);
    else if (high <= 0) result = positiveNormalTail(-high)-positiveNormalTail(-low);
    else result = 1-positiveNormalTail(high)-positiveNormalTail(-low);
    return Math.max(0, Math.min(1, result));
  }
  function readBounds(bounds, requireFinite) {
    if (!plainObject(bounds)) return null;
    const values = ['yMin','yMax','zMin','zMax'].map(key=>bounds[key]);
    if (values.some(v=>typeof v!=='number'||Number.isNaN(v)||(requireFinite&&!finite(v)))) return null;
    const [yMin,yMax,zMin,zMax] = values;
    if (yMin>yMax||zMin>zMax) return null;
    return {yMin,yMax,zMin,zMax};
  }
  function validateModel(model) {
    if (!plainObject(model)) return 'A density model object is required.';
    if (model.status !== 'ok' && model.status !== 'degenerate') return 'Density model is not usable; inspect its build status.';
    if (model.schema !== 'isotropic-gaussian-kde-1' || !Array.isArray(model.centers) ||
        !Array.isArray(model.weights) || model.weights.length!==model.centers.length ||
        !Array.isArray(model.issues)||model.issues.some(issue=>typeof issue!=='string')) return 'Malformed or unsupported KDE model.';
    if (model.centers.some(p=>!plainObject(p)||!finite(p.y)||!finite(p.z)) || model.weights.some(w=>!finite(w)||w<0)) return 'Invalid centers or weights.';
    let weightSum=0,compensation=0;
    for(const weight of model.weights){const corrected=weight-compensation,next=weightSum+corrected;compensation=(next-weightSum)-corrected;weightSum=next;}
    if (Math.abs(weightSum-(model.centers.length?1:0))>1e-12) return 'Conditional KDE weights must sum to one.';
    if (model.planeArrivalMass !== null) {
      const a = model.planeArrivalMass;
      if (!plainObject(a)||!Number.isSafeInteger(a.nTotal)||a.nTotal<1||a.nArrival!==model.centers.length||a.nArrival>a.nTotal||
          !finite(a.fraction)||Math.abs(a.fraction-a.nArrival/a.nTotal)>1e-15) return 'Invalid all-trial arrival accounting.';
    }
    if(model.status==='degenerate') return model.h===null?null:'Degenerate models must not claim a usable smoothing bandwidth.';
    if(model.centers.length<2||!finite(model.h)||model.h<=0||!finite(2*Math.PI*model.h*model.h)||!finite(1/(2*Math.PI*model.h*model.h))) return 'Malformed or unsupported KDE bandwidth.';
    if(model.supportZMin!==null&&!finite(model.supportZMin)) return 'supportZMin must be null or finite metres.';
    if(!Array.isArray(model.kernelNormalizers)||model.kernelNormalizers.length!==model.centers.length||model.kernelNormalizers.some(n=>!finite(n)||n<=0||n>1)) return 'Invalid per-kernel support normalizers.';
    return null;
  }
  function buildDensity(points, options = {}) {
    if (!Array.isArray(points)||!plainObject(options)) return failure('points must be an array and options an object.');
    if (points.some(p=>!plainObject(p)||!finite(p.y)||!finite(p.z))) return failure('Every point must contain finite y,z in metres; invalid samples are not silently dropped.');
    const allowed = ['sigma_h','minBandwidth','planeArrivalMass','supportZMin'];
    const unknown = Object.keys(options).filter(key=>!allowed.includes(key));
    if (unknown.length) return failure('Unknown options: '+unknown.join(', '));
    const minBandwidth = options.minBandwidth === undefined ? 0.001 : options.minBandwidth;
    const supportZMin = options.supportZMin === undefined ? null : options.supportZMin;
    if(supportZMin!==null&&!finite(supportZMin)) return failure('supportZMin must be null or a finite height in metres.');
    const manual = options.sigma_h !== undefined && options.sigma_h !== null;
    if (!finite(minBandwidth)||minBandwidth<=0 || (manual&&(!finite(options.sigma_h)||options.sigma_h<=0))) return failure('sigma_h and minBandwidth must be positive finite lengths in metres.');
    const centers = points.map(p=>({y:p.y,z:p.z})), n = centers.length;
    let arrival = null;
    if (options.planeArrivalMass !== undefined && options.planeArrivalMass !== null) {
      const a = options.planeArrivalMass;
      if (!plainObject(a)||!Number.isSafeInteger(a.nArrival)||!Number.isSafeInteger(a.nTotal)||
          a.nTotal<1||a.nArrival<0||a.nArrival!==n||a.nArrival>a.nTotal) return failure('planeArrivalMass needs nArrival=points.length and nTotal>=nArrival, with a positive all-trial integer denominator.');
      arrival = {nArrival:a.nArrival,nTotal:a.nTotal,fraction:a.nArrival/a.nTotal,
        interpretation:'Observed plane-arrival fraction among every attempted trial, including fit/integration failures. It is not true arrival or goal probability.'};
    }
    const result = {version:VERSION,schema:'isotropic-gaussian-kde-1',status:'degenerate',issues:[],sampleCount:n,
      centers,weights:n?centers.map(()=>1/n):[],h:null,supportZMin,kernelNormalizers:null,
      units:{coordinates:'m',h:'m',density:'m^-2',integratedMass:'dimensionless'},
      planeArrivalMass:arrival,
      bandwidth:{method:manual?'user_isotropic':'Scott_2D_isotropic',requested:manual?options.sigma_h:null,
        minBandwidth,raw:null,floorApplied:false,formula:'n^(-1/6) * sqrt((sampleVarianceY + sampleVarianceZ)/2)',
        interpretation:'Bandwidth is only a smoothing scale; it is not measurement noise or additional physical uncertainty.'},
      normalization:{support:supportZMin===null?'R^2':'R x [supportZMin, +Infinity)',conditionalIntegral:null,observedMassIntegral:null,truncatedToGoal:false,
        kernelRule:supportZMin===null?'Untruncated kernels':'Each kernel is truncated at supportZMin and separately renormalized; original mixture weights and bandwidth are unchanged.'},
      interpretation:{conditional:'Estimated density conditional on the supplied valid plane-arrival samples.',
        observedMass:'Conditional density times the observed arrival fraction; unavailable without all-trial counts.',
        warning:'Simulation/algorithm-conditional estimate, not calibrated real-world probability or a posterior.'}};
    if (n<2) {
      result.issues.push('At least two arrival samples are required. A zero/single-sample outcome is not an uncertainty density.');
      return result;
    }
    let meanY=0,meanZ=0,m2Y=0,m2Z=0;
    centers.forEach((p,i)=>{const dy=p.y-meanY,dz=p.z-meanZ;meanY+=dy/(i+1);meanZ+=dz/(i+1);m2Y+=dy*(p.y-meanY);m2Z+=dz*(p.z-meanZ);});
    const varianceY=m2Y/(n-1), varianceZ=m2Z/(n-1);
    result.sampleMoments={meanY,meanZ,varianceY,varianceZ,varianceDenominator:n-1};
    if (![meanY,meanZ,varianceY,varianceZ].every(finite) || varianceY<0||varianceZ<0) return failure('Sample moments exceed numerical range.');
    if (centers.every(p=>p.y===centers[0].y&&p.z===centers[0].z)) {
      result.issues.push('All samples coincide. Even a chosen smoothing width must not be displayed as inferred uncertainty.');
      return result;
    }
    const raw = manual ? options.sigma_h : Math.pow(n,-1/6)*Math.sqrt(varianceY/2+varianceZ/2);
    const h = Math.max(raw,minBandwidth);
    if (!finite(h)||h<=0||!finite(2*Math.PI*h*h)||!finite(1/(2*Math.PI*h*h))) return failure('Bandwidth is outside the representable numeric range.');
    const kernelNormalizers=centers.map(p=>supportZMin===null?1:normalInterval((supportZMin-p.z)/h,Infinity));
    if(kernelNormalizers.some(v=>v<=0||!finite(v)||!finite(1/(2*Math.PI*h*h*v)))) return failure('The retained kernel support mass is outside the representable numeric range.');
    result.status='ok';result.h=h;result.bandwidth.raw=raw;result.bandwidth.floorApplied=raw<minBandwidth;
    result.kernelNormalizers=kernelNormalizers;
    result.normalization.conditionalIntegral=1;
    result.normalization.observedMassIntegral=arrival===null?null:arrival.fraction;
    if (raw<minBandwidth) result.issues.push('The declared numerical bandwidth floor was applied; it adds smoothing, not evidence for uncertainty.');
    return result;
  }
  function densityValue(model,y,z) {
    if(model.supportZMin!==null&&z<model.supportZMin)return 0;
    let sum=0;
    for(let i=0;i<model.centers.length;i++) {
      const dy=(y-model.centers[i].y)/model.h,dz=(z-model.centers[i].z)/model.h;
      sum+=model.weights[i]*Math.exp(-0.5*(dy*dy+dz*dz))/model.kernelNormalizers[i];
    }
    return sum/(2*Math.PI*model.h*model.h);
  }
  function evaluateDensity(model,y,z) {
    const error=validateModel(model);if(error)return failure(error);
    if (!finite(y)||!finite(z)) return failure('Evaluation coordinates must be finite metres.');
    if(model.status==='degenerate')return {status:'degenerate',conditionalDensity:null,observedMassDensity:null,issues:model.issues.slice()};
    const value=densityValue(model,y,z);
    return {status:'ok',conditionalDensity:value,observedMassDensity:model.planeArrivalMass===null?null:value*model.planeArrivalMass.fraction,units:'m^-2'};
  }
  function integrateRectangle(model,bounds) {
    const error=validateModel(model);if(error)return failure(error);
    const b=readBounds(bounds,false);if(!b)return failure('Rectangle needs ordered numeric yMin,yMax,zMin,zMax; infinities are allowed.');
    if(model.status==='degenerate')return {status:'degenerate',conditionalMass:null,observedMass:null,issues:model.issues.slice()};
    let value=0;
    const zLower=model.supportZMin===null?b.zMin:Math.max(b.zMin,model.supportZMin);
    for(let i=0;i<model.centers.length;i++) {
      const p=model.centers[i];
      value+=model.weights[i]*normalInterval((b.yMin-p.y)/model.h,(b.yMax-p.y)/model.h)*normalInterval((zLower-p.z)/model.h,(b.zMax-p.z)/model.h)/model.kernelNormalizers[i];
    }
    value=Math.max(0,Math.min(1,value));
    return {status:'ok',bounds:b,conditionalMass:value,observedMass:model.planeArrivalMass===null?null:value*model.planeArrivalMass.fraction,
      units:'dimensionless',method:'Separable normal CDF intervals, respecting per-kernel supportZMin normalization; no renormalization to the queried rectangle.'};
  }
  function gridDensity(model,bounds,nx=61,ny=61) {
    const error=validateModel(model);if(error)return failure(error);
    const b=readBounds(bounds,true);
    if(!b||b.yMin===b.yMax||b.zMin===b.zMax||!finite(b.yMax-b.yMin)||!finite(b.zMax-b.zMin)) return failure('Grid requires finite, nonzero rectangle dimensions.');
    if(!Number.isSafeInteger(nx)||!Number.isSafeInteger(ny)||nx<2||ny<2||nx*ny>1000000) return failure('Grid dimensions must be integers >=2, with at most 1,000,000 points.');
    if(model.status==='degenerate')return {status:'degenerate',conditionalDensity:null,observedMassDensity:null,issues:model.issues.slice()};
    const y=Array.from({length:nx},(_,i)=>b.yMin+(b.yMax-b.yMin)*(i/(nx-1)));
    const z=Array.from({length:ny},(_,i)=>b.zMin+(b.zMax-b.zMin)*(i/(ny-1)));
    const conditionalDensity=z.map(zv=>y.map(yv=>densityValue(model,yv,zv)));
    return {status:'ok',bounds:b,nx,ny,y,z,layout:'density[zIndex][yIndex]',conditionalDensity,
      observedMassDensity:model.planeArrivalMass===null?null:conditionalDensity.map(row=>row.map(v=>v*model.planeArrivalMass.fraction)),
      rectangleMass:integrateRectangle(model,b),units:{coordinates:'m',density:'m^-2'},
      note:'Finite plotting window only; values are not renormalized to the grid or goal opening.'};
  }
  return {version:VERSION,buildDensity,evaluateDensity,integrateRectangle,gridDensity,normalCDF};
}
if(typeof module!=='undefined'&&module.exports)module.exports={createDensity};
if(typeof globalThis!=='undefined')globalThis.UPCDensity=createDensity();
