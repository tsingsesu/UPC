'use strict';
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const {spawnSync}=require('node:child_process');
const {createDensity}=require('../density.js');
const D=createDensity(),checks=[];
const check=(name,fn)=>{fn();checks.push({name,passed:true});};
const close=(actual,expected,tolerance=1e-12)=>assert(Math.abs(actual-expected)<=tolerance,`${actual} != ${expected} (tolerance ${tolerance})`);
const points=[{y:0,z:0},{y:.4,z:-.2},{y:1.1,z:.7}];
const model=D.buildDensity(points,{sigma_h:.35,planeArrivalMass:{nArrival:3,nTotal:5}});
const supported=D.buildDensity(points,{sigma_h:.35,supportZMin:.11,planeArrivalMass:{nArrival:3,nTotal:5}});
const all={yMin:-Infinity,yMax:Infinity,zMin:-Infinity,zMax:Infinity};

check('Automatic Scott bandwidth uses 2D sample variance and metre units',()=>{
 const p=[{y:-1,z:-2},{y:0,z:0},{y:1,z:2}];
 const m=D.buildDensity(p);assert.equal(m.status,'ok');
 close(m.h,Math.pow(3,-1/6)*Math.sqrt((1+4)/2));
 assert.equal(m.sampleMoments.varianceDenominator,2);assert.equal(m.units.density,'m^-2');
});
check('Numerical bandwidth floor is explicit and only smoothing',()=>{
 const m=D.buildDensity([{y:0,z:0},{y:1e-10,z:0}],{minBandwidth:.002});
 assert.equal(m.h,.002);assert.equal(m.bandwidth.floorApplied,true);assert(m.issues.length);
 const requested=D.buildDensity(points,{sigma_h:.01,minBandwidth:.02});
 assert.equal(requested.h,.02);assert.equal(requested.bandwidth.requested,.01);
});
check('Empty, one point, and repeated noiseless samples stay degenerate',()=>{
 for(const p of [[],[{y:1,z:2}],[{y:1,z:2},{y:1,z:2}]]){
  const m=D.buildDensity(p,{sigma_h:.1});assert.equal(m.status,'degenerate');
  assert.equal(m.h,null);assert.equal(D.evaluateDensity(m,1,2).conditionalDensity,null);
  assert.equal(D.integrateRectangle(m,all).conditionalMass,null);
 }
 const allFailed=D.buildDensity([],{planeArrivalMass:{nArrival:0,nTotal:10}});
 assert.equal(allFailed.status,'degenerate');assert.equal(allFailed.planeArrivalMass.fraction,0);
});
check('Invalid data and denominator accounting fail explicitly without silent filtering',()=>{
 const bads=[D.buildDensity([{y:NaN,z:1}]),D.buildDensity(points,{sigma_h:0}),D.buildDensity(points,{minBandwidth:-1}),
  D.buildDensity(points,{planeArrivalMass:{nArrival:2,nTotal:4}}),D.buildDensity(points,{planeArrivalMass:{nArrival:3,nTotal:2}}),
  D.buildDensity(points,{planeArrivalMass:{nArrival:3,nTotal:3.5}}),D.buildDensity(points,{sigma:1}),
  D.buildDensity(points,{sigma_h:1e-200,minBandwidth:1e-200}),D.buildDensity(points,{sigma_h:1e200})];
 bads.forEach(m=>assert.equal(m.status,'invalid_input'));
});
check('Full plane normalization and observed mass retain failed-trial denominator',()=>{
 assert.equal(model.status,'ok');close(model.weights.reduce((a,b)=>a+b),1);
 const integral=D.integrateRectangle(model,all);close(integral.conditionalMass,1);close(integral.observedMass,.6);
 const value=D.evaluateDensity(model,.2,.1);close(value.observedMassDensity,.6*value.conditionalDensity);
 const unaccounted=D.buildDensity(points,{sigma_h:.35});
 assert.equal(D.evaluateDensity(unaccounted,0,0).observedMassDensity,null);
 assert.equal(D.integrateRectangle(unaccounted,all).observedMass,null);
});
check('Explicit mixture PDF agrees with direct formula and translation invariance',()=>{
 const y=.2,z=.1,expected=points.reduce((s,p)=>s+Math.exp(-((y-p.y)**2+(z-p.z)**2)/(2*.35**2))/(2*Math.PI*.35**2),0)/3;
 close(D.evaluateDensity(model,y,z).conditionalDensity,expected);
 const moved=D.buildDensity(points.map(p=>({y:p.y+10,z:p.z-4})),{sigma_h:.35});
 close(D.evaluateDensity(moved,y+10,z-4).conditionalDensity,expected);
});
check('Density and bandwidth obey metre-to-scaled-length dimensional transformation',()=>{
 const scale=100,scaled=D.buildDensity(points.map(p=>({y:p.y*scale,z:p.z*scale})),{sigma_h:model.h*scale});
 close(D.evaluateDensity(scaled,20,10).conditionalDensity*scale**2,D.evaluateDensity(model,.2,.1).conditionalDensity);
});
check('Rectangle mass is not renormalized inside a goal opening or display window',()=>{
 const b={yMin:-3.66,yMax:3.66,zMin:0,zMax:2.44},m=D.integrateRectangle(model,b);
 assert(m.conditionalMass>0&&m.conditionalMass<.9);close(m.observedMass,.6*m.conditionalMass);
 close(D.integrateRectangle(model,{yMin:0,yMax:0,zMin:0,zMax:1}).conditionalMass,0);
 const left=D.integrateRectangle(model,{yMin:-Infinity,yMax:.2,zMin:-Infinity,zMax:Infinity});
 const right=D.integrateRectangle(model,{yMin:.2,yMax:Infinity,zMin:-Infinity,zMax:Infinity});
 close(left.conditionalMass+right.conditionalMass,1);
});
check('Grid has shared coordinates, z-row/y-column layout, and no local normalization',()=>{
 const b={yMin:-1,yMax:2,zMin:-2,zMax:2},g=D.gridDensity(model,b,31,41);
 assert.equal(g.status,'ok');assert.equal(g.conditionalDensity.length,41);assert.equal(g.conditionalDensity[0].length,31);
 close(g.conditionalDensity[20][15],D.evaluateDensity(model,g.y[15],g.z[20]).conditionalDensity);
 close(g.observedMassDensity[20][15],g.conditionalDensity[20][15]*.6);
 assert(g.rectangleMass.conditionalMass<1);assert.equal(D.gridDensity(model,b,1,4).status,'invalid_input');
 assert.equal(D.gridDensity(model,all,3,3).status,'invalid_input');
});
check('Wide-grid quadrature approaches the global normalized density',()=>{
 const g=D.gridDensity(model,{yMin:-4,yMax:5,zMin:-4,zMax:5},301,301);
 let sum=0;for(let j=0;j<g.ny;j++)for(let i=0;i<g.nx;i++)sum+=g.conditionalDensity[j][i]*(i===0||i===g.nx-1?.5:1)*(j===0||j===g.ny-1?.5:1);
 sum*=(g.y[1]-g.y[0])*(g.z[1]-g.z[0]);close(sum,1,2e-12);
});
check('JSON roundtrip preserves exported function, complete centers, weights and semantics',()=>{
 const restored=JSON.parse(JSON.stringify(model));assert.deepEqual(D.evaluateDensity(restored,.2,.1),D.evaluateDensity(model,.2,.1));
 assert.deepEqual(D.integrateRectangle(restored,all),D.integrateRectangle(model,all));
 const damaged=JSON.parse(JSON.stringify(model));damaged.weights[0]=0;assert.equal(D.evaluateDensity(damaged,0,0).status,'invalid_input');
 assert.equal(D.evaluateDensity({status:'degenerate'},0,0).status,'invalid_input');
 const malformedWidth=JSON.parse(JSON.stringify(model));malformedWidth.h=1e200;assert.equal(D.evaluateDensity(malformedWidth,0,0).status,'invalid_input');
 const input=JSON.stringify(points),options={sigma_h:.35};D.buildDensity(points,options);assert.equal(JSON.stringify(points),input);assert.deepEqual(options,{sigma_h:.35});
});
check('Factory serializes into an isolated Worker-like context without dependencies',()=>{
 const W=vm.runInNewContext('('+createDensity.toString()+')()');
 close(W.evaluateDensity(JSON.parse(JSON.stringify(model)),.2,.1).conditionalDensity,D.evaluateDensity(model,.2,.1).conditionalDensity);
 const exported='const model='+JSON.stringify(model)+'; const api=('+createDensity.toString()+')(); const f=(y,z)=>api.evaluateDensity(model,y,z).conditionalDensity; f(.2,.1);';
 close(vm.runInNewContext(exported),D.evaluateDensity(model,.2,.1).conditionalDensity);
});
check('Physical lower support truncates and normalizes every kernel without changing h or y support',()=>{
 assert.equal(supported.status,'ok');assert.equal(supported.h,model.h);assert.deepEqual(supported.weights,model.weights);
 close(D.integrateRectangle(supported,all).conditionalMass,1);
 close(D.integrateRectangle(supported,{...all,zMin:.11}).conditionalMass,1);
 close(D.integrateRectangle(supported,{...all,zMax:.11}).conditionalMass,0);
 close(D.integrateRectangle(supported,all).observedMass,.6);
 assert.equal(D.evaluateDensity(supported,0,.109999).conditionalDensity,0);
 assert(D.evaluateDensity(supported,0,.11).conditionalDensity>0);
 assert(D.evaluateDensity(supported,5,.11).conditionalDensity>0);
 assert(supported.kernelNormalizers.every(v=>v>0&&v<1));
 assert.equal(D.buildDensity(points,{supportZMin:Infinity}).status,'invalid_input');
});
check('Truncated function JSON and standalone manual formula retain every normalizer and boundary',()=>{
 const m=JSON.parse(JSON.stringify(supported));
 const direct=(y,z)=>z<m.supportZMin?0:m.centers.reduce((sum,p,i)=>sum+m.weights[i]*Math.exp(-((y-p.y)**2+(z-p.z)**2)/(2*m.h*m.h))/m.kernelNormalizers[i],0)/(2*Math.PI*m.h*m.h);
 for(const [y,z] of [[.2,.1],[0,.11],[1,.5],[3,2]])close(direct(y,z),D.evaluateDensity(m,y,z).conditionalDensity);
 const W=vm.runInNewContext('('+createDensity.toString()+')()');close(W.evaluateDensity(m,.2,.2).conditionalDensity,direct(.2,.2));
});

// Independent numerical reference: SciPy's normal CDF/PDF and adaptive
// double quadrature, not a transcription of the JS CDF approximation.
const referenceRequest={points,h:.35,cdfInputs:[-30,-10,-8,-5,-3,-2.99,-1,0,1,2.99,3,5,8,10,30],
 supportZMin:.11,supportEvaluationPoints:[[0,.1],[0,.11],[.2,.2],[1,1]],
 supportRectangles:[[-1,1,-1,1],[-3.66,3.66,0,2.44],[-1,1,-1,.1],[-1,1,.11,.5]],
 evaluationPoints:[[0,0],[.2,.1],[4,-2],[-8,0]],
 rectangles:[[-1,1,-1,1],[-3.66,3.66,0,2.44],[3,3.1,-.2,.4],[8,8.1,-.2,.4],[-8,-7.9,-.2,.4],[.01,.010001,-.01,.010001]]};
const python=process.env.SCIPY_PYTHON||'python';
// Python must provide NumPy/SciPy; set SCIPY_PYTHON to select its executable.
const pythonScript=String.raw`
import json, sys, math, scipy
from scipy.stats import norm
from scipy.integrate import dblquad
r=json.load(sys.stdin)
p=r['points']; h=r['h']
def interval(lo,hi):
    if lo>=0: return float(norm.sf(lo)-norm.sf(hi))
    return float(norm.cdf(hi)-norm.cdf(lo))
def pdf(y,z):
    return sum(float(norm.pdf((y-a['y'])/h)*norm.pdf((z-a['z'])/h)/h**2) for a in p)/len(p)
masses=[]
for yl,yu,zl,zu in r['rectangles']:
    masses.append(sum(interval((yl-a['y'])/h,(yu-a['y'])/h)*interval((zl-a['z'])/h,(zu-a['z'])/h) for a in p)/len(p))
q,err=dblquad(lambda z,y:pdf(y,z),-1,1,lambda y:-1,lambda y:1,epsabs=1e-11,epsrel=1e-11)
zmin=r['supportZMin']; normalizers=[float(norm.sf((zmin-a['z'])/h)) for a in p]
def support_pdf(y,z):
    if z<zmin: return 0.0
    return sum(float(norm.pdf((y-a['y'])/h)*norm.pdf((z-a['z'])/h)/h**2)/normalizers[i] for i,a in enumerate(p))/len(p)
support_masses=[]
for yl,yu,zl,zu in r['supportRectangles']:
    zl=max(zl,zmin)
    support_masses.append(0.0 if zl>=zu else sum(interval((yl-a['y'])/h,(yu-a['y'])/h)*interval((zl-a['z'])/h,(zu-a['z'])/h)/normalizers[i] for i,a in enumerate(p))/len(p))
support_q,support_err=dblquad(lambda z,y:support_pdf(y,z),-5,6,lambda y:zmin,lambda y:6,epsabs=1e-10,epsrel=1e-10)
print(json.dumps(dict(scipyVersion=scipy.__version__,method='scipy.stats.norm cdf/sf/pdf + scipy.integrate.dblquad',cdf=[float(norm.cdf(x)) for x in r['cdfInputs']],density=[pdf(y,z) for y,z in r['evaluationPoints']],masses=masses,quadrature=q,quadratureErrorEstimate=err,supportNormalizers=normalizers,supportDensity=[support_pdf(y,z) for y,z in r['supportEvaluationPoints']],supportMasses=support_masses,supportQuadrature=support_q,supportQuadratureErrorEstimate=support_err)))
`;
const processResult=spawnSync(python,['-c',pythonScript],{input:JSON.stringify(referenceRequest),encoding:'utf8',windowsHide:true});
assert.equal(processResult.status,0,processResult.stderr);const reference=JSON.parse(processResult.stdout);
let maxCdfAbsError=0,maxDensityAbsError=0,maxRectangleAbsError=0,maxTailRelativeError=0;
check('Normal CDF independently agrees with SciPy across central and far-tail inputs',()=>{
 referenceRequest.cdfInputs.forEach((x,i)=>{const e=Math.abs(D.normalCDF(x)-reference.cdf[i]);maxCdfAbsError=Math.max(maxCdfAbsError,e);assert(e<8e-15);});
 close(D.normalCDF(-Infinity),0);close(D.normalCDF(Infinity),1);
});
check('KDE PDF independently agrees with SciPy normal products',()=>{
 referenceRequest.evaluationPoints.forEach(([y,z],i)=>{const e=Math.abs(D.evaluateDensity(model,y,z).conditionalDensity-reference.density[i]);maxDensityAbsError=Math.max(maxDensityAbsError,e);assert(e<1e-14);});
});
check('Rectangles agree with independent SciPy CDF and adaptive double quadrature, including tiny tails',()=>{
 referenceRequest.rectangles.forEach(([yMin,yMax,zMin,zMax],i)=>{
  const actual=D.integrateRectangle(model,{yMin,yMax,zMin,zMax}).conditionalMass,expected=reference.masses[i],e=Math.abs(actual-expected);
  maxRectangleAbsError=Math.max(maxRectangleAbsError,e);assert(e<2e-14);
  if(i===2||i===3||i===4){const relative=e/expected;maxTailRelativeError=Math.max(maxTailRelativeError,relative);assert(relative<1e-10);assert(actual>0);}
 });
 close(D.integrateRectangle(model,{yMin:-1,yMax:1,zMin:-1,zMax:1}).conditionalMass,reference.quadrature,2e-12);
});
check('Truncated kernels independently agree with SciPy at boundary and in rectangle integrals',()=>{
 reference.supportNormalizers.forEach((value,i)=>close(supported.kernelNormalizers[i],value,2e-15));
 referenceRequest.supportEvaluationPoints.forEach(([y,z],i)=>close(D.evaluateDensity(supported,y,z).conditionalDensity,reference.supportDensity[i],2e-14));
 referenceRequest.supportRectangles.forEach(([yMin,yMax,zMin,zMax],i)=>close(D.integrateRectangle(supported,{yMin,yMax,zMin,zMax}).conditionalMass,reference.supportMasses[i],2e-14));
 close(reference.supportQuadrature,1,2e-10);
 close(D.integrateRectangle(supported,{yMin:-5,yMax:6,zMin:.11,zMax:6}).conditionalMass,reference.supportQuadrature,2e-10);
});
const report={passed:true,checks:checks.length,version:D.version,node:process.version,tests:checks,
 independentReference:{...reference,request:referenceRequest},
 errors:{maxCdfAbsError,maxDensityAbsError,maxRectangleAbsError,maxTailRelativeError},
 interpretation:'Numerical KDE verification only; no model ranking, physical validation, or calibrated uncertainty coverage.'};
if(process.argv.includes('--write-report'))fs.writeFileSync(path.join(__dirname,'density.validation.json'),JSON.stringify(report,null,2));
console.log(JSON.stringify({passed:true,checks:report.checks,scipy:reference.scipyVersion,errors:report.errors},null,2));
