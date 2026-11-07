/* Fixed-design synthetic robustness audit. Run with --write-report to save JSON/CSV.
   The truth generator uses the SAME RK4 implementation at finer resolution.
   It is not an independent physical experiment or an independent algorithm. */
'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const P=require('./physics.js'),started=Date.now();
const SEED=2026110717,SIGMA=.01,TRUTH_H=.0001,TRUTH_OBS_H=.00005;
function xorshift32(seed){let x=seed>>>0;return ()=>{x^=x<<13;x^=x>>>17;x^=x<<5;return (x>>>0)/4294967296;};}
function gaussian(r){return Math.sqrt(-2*Math.log(Math.max(r(),1e-15)))*Math.cos(2*Math.PI*r());}
const uniform=(r,a,b)=>a+(b-a)*r();
function wilson(k,n){if(!n)return null;const z=1.959963984540054,p=k/n,d=1+z*z/n,center=(p+z*z/(2*n))/d,radius=z*Math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;
  return {successes:k,denominator:n,frequency:p,standardError:Math.sqrt(p*(1-p)/n),nominalConfidenceLevel:.95,wilsonInterval:[Math.max(0,center-radius),Math.min(1,center+radius)],
    meaning:'有限独立重复样本下，对当前生成算法稳定性频率的Wilson区间；不是现实进球概率区间，也不是落点95%覆盖区域'};}
function groupSummary(rows){const statuses={},truthStatuses={};for(const r of rows){statuses[r.estimateStatus]=(statuses[r.estimateStatus]||0)+1;truthStatuses[r.truthStatus]=(truthStatuses[r.truthStatus]||0)+1;}
  const paired=rows.filter(r=>r.errors!==null),n=paired.length;
  const mean=key=>n?paired.reduce((s,r)=>s+r.errors[key],0)/n:null;
  const rms=key=>n?Math.sqrt(paired.reduce((s,r)=>s+r.errors[key]**2,0)/n):null;
  return {total:rows.length,estimateStatusCounts:statuses,truthStatusCounts:truthStatuses,fitFailed:rows.filter(r=>!r.fitOk).length,
    boundaryHits:rows.filter(r=>r.boundaryHit).length,residualWarnings:rows.filter(r=>r.compatibility&&r.compatibility.flag!=='not_flagged_by_residual_rule').length,
    pairedGoalErrorDenominator:n,excludedFromGoalErrorButRetained:rows.length-n,
    meanError_y_m:mean('y_m'),meanError_z_m:mean('z_m'),RMSE_y_m:rms('y_m'),RMSE_z_m:rms('z_m'),RMSE_plane_m:n?Math.sqrt(paired.reduce((s,r)=>s+r.errors.y_m**2+r.errors.z_m**2,0)/n):null,
    RMSE_time_s:rms('time_s'),meanTimeError_s:mean('time_s'),
    note:'落点与时间误差只对真值和预测都到达球门的配对样本计算；全部失败仍保留在总数和样本表中。'};
}
function auditGroup(name,count,seed,mismatch){const rows=[];
  for(let i=0;i<count;i++){const sampleSeed=(seed+Math.imul(i+1,0x9e3779b1))>>>0,r=xorshift32(sampleSeed);
    const th=[uniform(r,24,36),uniform(r,2,8),uniform(r,4,9),uniform(r,.15,.45),uniform(r,-.25,.25),uniform(r,-.55,.15)];
    const wind=mismatch?[uniform(r,-5,5),(i%2?-1:1)*uniform(r,3,8),0]:[0,0,0];
    const truthConfig={...P.defaultConfig(),mode:'forward',velocity:th.slice(0,3),Cd:th[3],b:[0,...th.slice(4)],wind,h:TRUTH_H,fitH:TRUTH_OBS_H};
    const clean=P.predictPositions(truthConfig,th,[.05,.1]),truth=P.integrateTrajectory(truthConfig,th,{collect:false}).prediction;
    const noise=Array.from({length:6},()=>SIGMA*gaussian(r));
    const obs=clean.map((p,j)=>({t:[.05,.1][j],position:p.map((x,k)=>x+noise[3*j+k])}));
    const fitConfig={...P.defaultConfig(),mode:'inverse',sigma:SIGMA,observations:obs,wind:[0,0,0]};
    const estimate=P.solve(fitConfig),paired=truth.status==='crossed_goal_plane'&&estimate.prediction.status==='crossed_goal_plane';
    let noiselessMismatch=null;
    if(mismatch){const cleanFit=P.solve({...fitConfig,observations:clean.map((position,j)=>({t:[.05,.1][j],position}))});
      noiselessMismatch={prediction:cleanFit.prediction,fitRMSE_m:cleanFit.fit.rmse,boundaryHit:cleanFit.fit.boundaryHit,
        errors:truth.status==='crossed_goal_plane'&&cleanFit.prediction.status==='crossed_goal_plane'?{y_m:cleanFit.prediction.yGoal-truth.yGoal,z_m:cleanFit.prediction.zGoal-truth.zGoal,time_s:cleanFit.prediction.tGoal-truth.tGoal}:null};}
    rows.push({id:name+'_'+String(i+1).padStart(3,'0'),group:name,index:i+1,seed:sampleSeed,trueTheta:th,trueWind_m_s:wind,
      cleanObservations:clean,noise_m:noise,observations:obs,trueConfig:truthConfig,fitConfig,truthStatus:truth.status,truth,
      estimateStatus:estimate.prediction.status,estimate:estimate.prediction,estimatedTheta:estimate.theta,fitOk:estimate.fit.ok,
      fitRMSE_m:estimate.fit.rmse,fitResiduals_m:estimate.fit.residuals,boundaryHit:estimate.fit.boundaryHit,compatibility:estimate.fit.diagnostics.compatibility||null,
      optimizerStarts:estimate.fit.starts.map(s=>({converged:s.converged,accepted:s.accepted,iterations:s.iterations,reason:s.reason,boundaryHit:s.boundaryHit,sse:s.sse})),
      errors:paired?{y_m:estimate.prediction.yGoal-truth.yGoal,z_m:estimate.prediction.zGoal-truth.zGoal,time_s:estimate.prediction.tGoal-truth.tGoal}:null,noiselessMismatch});
  }
  return {name,sampleCount:count,seed,summary:groupSummary(rows),samples:rows};
}
const sameModel=auditGroup('same_model_M3',60,SEED,false);
const mismatch=auditGroup('omitted_wind_M3',30,SEED+1,true);
const checks=[];
function check(name,fn){try{const evidence=fn();checks.push({name,passed:true,evidence:evidence||null});}catch(e){checks.push({name,passed:false,error:e.stack});}}
check('all 90 predefined samples retained with reconstructable inputs',()=>{
  assert.equal(sameModel.samples.length,60);assert.equal(mismatch.samples.length,30);
  for(const g of [sameModel,mismatch]){assert.equal(Object.values(g.summary.estimateStatusCounts).reduce((a,b)=>a+b,0),g.sampleCount);assert.equal(g.summary.pairedGoalErrorDenominator+g.summary.excludedFromGoalErrorButRetained,g.sampleCount);
    for(const s of g.samples){assert.equal(s.noise_m.length,6);assert.equal(s.observations.length,2);assert.ok(s.estimateStatus);if(s.estimateStatus!=='crossed_goal_plane')assert.equal(s.errors,null);}}
});
check('near-tangent directed event exposes conditioning',()=>{
  const c={...P.defaultConfig(),model:'M1',mode:'forward',g:0,initial:[0,0,2],velocity:[2,0,0],constantAcceleration:[-2,0,0],goalX:.99999999};
  const r=P.solve(c),exact=1-Math.sqrt(1-c.goalX);assert.equal(r.prediction.status,'crossed_goal_plane');assert.equal(r.prediction.nearTangent,true);assert.ok(Math.abs(r.prediction.tGoal-exact)<1e-7);
  return {prediction:r.prediction,analyticTime_s:exact,error_s:Math.abs(r.prediction.tGoal-exact),note:'接近切线时根的条件性差；只测试此指定例，不把普通误差阈值推广到所有切线工况。'};
});
check('legal numerical extremes return finite trajectories or explicit states',()=>{
  const base={...P.defaultConfig(),mode:'forward'};
  const cases=[{id:'fast_kick',config:{...base,velocity:[90,0,30],Cd:.6,b:[0,1.5,-1.5]}},
    {id:'strong_effective_b',config:{...base,initial:[0,0,2],velocity:[50,0,5],b:[0,7,7]}},
    {id:'strong_headwind',config:{...base,wind:[-50,0,0],velocity:[30,0,6]}},
    {id:'weak_launch_ground',config:{...base,velocity:[1,0,.1]}},
    {id:'no_crossing',config:{...base,model:'M0',g:0,initial:[0,0,2],velocity:[0,0,0]}}];
  return cases.map(({id,config})=>{const r=P.solve(config);assert.equal(r.issues.length,0);assert.ok(['crossed_goal_plane','ground_before_goal','no_crossing_within_horizon','integration_failed','outside_model_domain'].includes(r.prediction.status));
    assert.ok(r.trajectory.every(q=>[q.t,...q.position,...q.velocity,...q.acceleration].every(Number.isFinite)));return {id,prediction:r.prediction,pointCount:r.trajectory.length};});
});
check('incompatible observations preserve residuals and trigger explicit warning',()=>{
  const c=P.defaultConfig();c.observations[0].position[0]=20;c.observations[1].position[0]=40;
  const r=P.solve(c);assert.equal(r.fit.diagnostics.compatibility.flag,'incompatible_with_declared_noise');assert.ok(r.issues.length);assert.ok(r.fit.boundaryHit);
  return {fitOkMeansOptimizerSuccess:r.fit.ok,rmse_m:r.fit.rmse,compatibility:r.fit.diagnostics.compatibility,issues:r.issues};
});
check('near-coincident frames report severe conditioning or failure',()=>{
  const t1=.05,t2=.05000001,th=[30,5,7,.3,.1,-.2],c={...P.defaultConfig(),mode:'forward',velocity:th.slice(0,3),Cd:th[3],b:[0,...th.slice(4)],fitH:TRUTH_OBS_H};
  const positions=P.predictPositions(c,th,[t1,t2]);const r=P.solve({...P.defaultConfig(),observations:positions.map((position,i)=>({t:[t1,t2][i],position}))});
  if(r.fit.ok)assert.ok(r.fit.diagnostics.rank<6||r.fit.diagnostics.condition>1e6);
  return {status:r.prediction.status,fitOk:r.fit.ok,rank:r.fit.diagnostics.rank||null,condition:r.fit.diagnostics.condition||null,note:'近同时观测缺少独立动力学信息，不要求唯一恢复真实参数。'};
});
const boot=P.uncertainty(P.defaultConfig(),200,2026110788,null,'bootstrap');
const bootstrapStability=[50,100,200].map(n=>{const sub=boot.samples.slice(0,n),arrivals=sub.filter(s=>s.status==='crossed_goal_plane');return {sampleCount:n,
  arrivalFrequency:wilson(arrivals.length,n),centerInsideAmongArrivals:wilson(arrivals.filter(s=>s.centerInside).length,arrivals.length),
  clearanceInsideAmongArrivals:wilson(arrivals.filter(s=>s.clearanceInside).length,arrivals.length),
  boundaryHits:sub.filter(s=>s.boundaryHit).length,fitFailed:sub.filter(s=>s.status==='fit_failed').length,
  warning:'三个前缀共享样本而相互相关；该对照不证明收敛，也不评估现实覆盖率。'};});
check('Wilson endpoints and conditional denominator accounting',()=>{
  assert.equal(wilson(0,0),null);assert.equal(wilson(0,10).wilsonInterval[0],0);assert.ok(wilson(10,10).wilsonInterval[1]<=1);
  for(const x of bootstrapStability){assert.equal(x.arrivalFrequency.denominator,x.sampleCount);assert.equal(x.centerInsideAmongArrivals.denominator,x.arrivalFrequency.successes);assert.ok(x.centerInsideAmongArrivals.wilsonInterval.every(Number.isFinite));}
});
const mismatchCleanRows=mismatch.samples.map(s=>({...s,estimateStatus:s.noiselessMismatch.prediction.status,errors:s.noiselessMismatch.errors}));
mismatch.noiselessStructuralError={RMSE_plane_m:groupSummary(mismatchCleanRows).RMSE_plane_m,denominator:mismatchCleanRows.filter(s=>s.errors).length,
  meaning:'相同30个带风真值的无噪声静风反演对照，分离测量噪声和所选模型失配；未代表所有可能失配机制。'};
const report={schemaVersion:'1.0',schemaRef:'浏览器鲁棒性验证.schema.json',artifactKind:'synthetic_robustness_validation',codeVersion:P.VERSION,seed:SEED,
  csvNullEncoding:'Empty cell means unavailable event/failed fit, never a measured zero.',
  provenance:{node:process.version,platform:process.platform,corePath:'src/web/physics.js',coreSHA256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'physics.js'))).digest('hex'),
    generator:'Same JavaScript classical RK4 implementation; smaller step sizes only, NOT algorithmically independent.',
    independentImplementationReferences:['data/验证基准/浏览器对照.json (MATLAB)', 'data/验证基准/浏览器内核测试.json (MATLAB comparison tests)'],
    realExperimentalValidation:false,newCoverageClaim:false},
  design:{sameModelCount:60,mismatchCount:30,initialPositionExact:true,observationTimes_s:[.05,.1],sigmaPerCoordinate_m:SIGMA,observationNoise:'independent zero-mean Gaussian; explicit hypothetical noise model',
    rng:'xorshift32 for truth parameter/noise design; Box–Muller; per-index stream seed retained',
    truthParameterDistribution:{vx_m_s:[24,36],vy_m_s:[2,8],vz_m_s:[4,9],Cd:[.15,.45],by_s_inv:[-.25,.25],bz_s_inv:[-.55,.15],distribution:'independent uniform for each listed interval; bx=0'},
    mismatch:'Truth wind x uniform [-5,5] m/s, |y| uniform [3,8] m/s with alternating sign, z=0; fit fixes all wind components to 0.',
    truthPredictionStep_s:TRUTH_H,truthObservationStep_s:TRUTH_OBS_H,fitObservationStep_s:.005,fitPredictionStep_s:.001,
    populations:'Two predeclared synthetic design populations, kept separate. Neither is a measured real-world football population.'},
  sampleSchema:{id:'string',group:'string',seed:'uint32',trueTheta:'[vx,vy,vz,Cd,by,bz] in SI',trueWind_m_s:'3-vector',cleanObservations:'2 by 3 m',noise_m:'[x1,y1,z1,x2,y2,z2] m',observations:'two {t,position} rows',trueConfig:'complete generator config',fitConfig:'complete fit config',truth:'prediction event object',estimate:'prediction event object',estimatedTheta:'six numbers or null',fitOk:'boolean optimizer success only',boundaryHit:'boolean overlapping other statuses',fitResiduals_m:'six residuals or empty on failure',compatibility:'engineering residual warning, not a formal statistical test',errors:'{y_m,z_m,time_s} only if truth and prediction both reach; otherwise null'},
  sameModel,mismatch,bootstrap:{sampleCount:200,seed:2026110788,method:boot.method,noiseAssumption:boot.noiseAssumption,summary:boot.summary,finiteSampleFrequencyAudit:bootstrapStability,samples:boot.samples},
  checks,passed:checks.every(c=>c.passed),elapsedMs:Date.now()-started,
  limitations:['No posterior probability or real goal probability is inferred.','No 95% prediction coverage is asserted.','Synthetic truth shares the RK4 code with the estimator; finer resolution cannot detect all common implementation errors.','Only omitted uniform wind is studied as model mismatch; other misspecifications can behave differently.','Boundary hits and all failure states remain in the denominators.']};
if(process.argv.includes('--write-report')){
  const dest=path.join(__dirname,'../../data/验证基准');fs.mkdirSync(dest,{recursive:true});fs.writeFileSync(path.join(dest,'浏览器鲁棒性验证.json'),JSON.stringify(report,null,2),'utf8');
  const header=['id','group','seed','truth_status','estimate_status','fit_ok','boundary_hit','fit_rmse_m','true_vx_m_s','true_vy_m_s','true_vz_m_s','true_Cd','true_by_s_inv','true_bz_s_inv','wind_x_m_s','wind_y_m_s','wind_z_m_s','obs_x1_m','obs_y1_m','obs_z1_m','obs_x2_m','obs_y2_m','obs_z2_m','estimated_vx_m_s','estimated_vy_m_s','estimated_vz_m_s','estimated_Cd','estimated_by_s_inv','estimated_bz_s_inv','true_yGoal_m','true_zGoal_m','true_tGoal_s','estimated_yGoal_m','estimated_zGoal_m','estimated_tGoal_s','error_y_m','error_z_m','error_t_s','compatibility_flag'];
  const q=x=>x===null||x===undefined?'':typeof x==='string'?'"'+x.replaceAll('"','""')+'"':String(x);
  const lines=[header.join(',')];for(const r of [...sameModel.samples,...mismatch.samples])lines.push([r.id,r.group,r.seed,r.truthStatus,r.estimateStatus,r.fitOk,r.boundaryHit,r.fitRMSE_m,...r.trueTheta,...r.trueWind_m_s,...r.observations.flatMap(o=>o.position),...(r.estimatedTheta||Array(6).fill(null)),r.truth.yGoal,r.truth.zGoal,r.truth.tGoal,r.estimate.yGoal,r.estimate.zGoal,r.estimate.tGoal,r.errors?.y_m,r.errors?.z_m,r.errors?.time_s,r.compatibility?.flag].map(q).join(','));
  fs.writeFileSync(path.join(dest,'浏览器鲁棒性样本.csv'),lines.join('\r\n')+'\r\n','utf8');
  const bootstrapHeader=['index','seed','status','boundary_hit','fit_rmse_m','yGoal_m','zGoal_m','tGoal_s','center_inside','clearance_inside','vx_m_s','vy_m_s','vz_m_s','Cd','by_s_inv','bz_s_inv','obs_x1_m','obs_y1_m','obs_z1_m','obs_x2_m','obs_y2_m','obs_z2_m'];
  const bootstrapLines=[bootstrapHeader.join(',')];for(const s of boot.samples)bootstrapLines.push([s.index,s.seed,s.status,s.boundaryHit,s.rmse,s.yGoal,s.zGoal,s.tGoal,s.centerInside,s.clearanceInside,...(s.theta||Array(6).fill(null)),...(s.observations?s.observations.flatMap(o=>o.position):Array(6).fill(null))].map(q).join(','));
  fs.writeFileSync(path.join(dest,'浏览器Bootstrap样本.csv'),bootstrapLines.join('\r\n')+'\r\n','utf8');
  const schema={$schema:'https://json-schema.org/draft/2020-12/schema',title:'UPC synthetic robustness report 1.0',type:'object',required:['schemaVersion','artifactKind','provenance','design','sameModel','mismatch','bootstrap','checks','passed'],
    properties:{schemaVersion:{const:'1.0'},artifactKind:{const:'synthetic_robustness_validation'},provenance:{type:'object',required:['coreSHA256','generator','realExperimentalValidation'],properties:{coreSHA256:{type:'string',pattern:'^[a-f0-9]{64}$'},realExperimentalValidation:{const:false}}},
      sameModel:{allOf:[{$ref:'#/$defs/group'},{properties:{sampleCount:{const:60},samples:{minItems:60,maxItems:60}}}]},mismatch:{allOf:[{$ref:'#/$defs/group'},{properties:{sampleCount:{const:30},samples:{minItems:30,maxItems:30}}}]},
      bootstrap:{type:'object',required:['sampleCount','samples','summary','finiteSampleFrequencyAudit'],properties:{sampleCount:{const:200},samples:{type:'array',minItems:200,maxItems:200}}},checks:{type:'array',items:{type:'object',required:['name','passed'],properties:{name:{type:'string'},passed:{type:'boolean'}}}},passed:{type:'boolean'}},
    $defs:{group:{type:'object',required:['name','sampleCount','seed','summary','samples'],properties:{name:{type:'string'},sampleCount:{type:'integer',minimum:1},seed:{type:'integer'},summary:{type:'object'},samples:{type:'array',items:{$ref:'#/$defs/sample'}}}},
      sample:{type:'object',required:['id','group','seed','trueTheta','trueWind_m_s','observations','noise_m','truthStatus','estimateStatus','fitOk','boundaryHit','errors'],properties:{id:{type:'string'},group:{type:'string'},seed:{type:'integer',minimum:0,maximum:4294967295},trueTheta:{type:'array',minItems:6,maxItems:6,items:{type:'number'}},trueWind_m_s:{type:'array',minItems:3,maxItems:3,items:{type:'number'}},noise_m:{type:'array',minItems:6,maxItems:6,items:{type:'number'}},observations:{type:'array',minItems:2,maxItems:2},truthStatus:{type:'string'},estimateStatus:{type:'string'},fitOk:{type:'boolean'},boundaryHit:{type:'boolean'},errors:{anyOf:[{type:'null'},{type:'object',required:['y_m','z_m','time_s'],properties:{y_m:{type:'number'},z_m:{type:'number'},time_s:{type:'number'}}}]}}}}};
  fs.writeFileSync(path.join(dest,'浏览器鲁棒性验证.schema.json'),JSON.stringify(schema,null,2),'utf8');
}
console.log(JSON.stringify({passed:report.passed,sameModel:report.sameModel.summary,mismatch:report.mismatch.summary,noiselessStructuralError:mismatch.noiselessStructuralError,bootstrapStability,checks,elapsedMs:report.elapsedMs},null,2));
if(!report.passed)process.exitCode=1;
