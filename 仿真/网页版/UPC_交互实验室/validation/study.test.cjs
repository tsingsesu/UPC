'use strict';
// Independent contract tests for paired perturbations, failure accounting and
// stability decisions. No model is assumed to be the physical ground truth.
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const vm=require('node:vm');
const {createDensity}=require('../density.js');
const {createStudy}=require('../study.js');
const physicsPath='../physics.js';
const {createPhysics}=require(physicsPath);
const P=createPhysics(),D=createDensity(),S=createStudy(P,D),checks=[];
const clone=x=>JSON.parse(JSON.stringify(x));
const close=(a,b,tolerance=1e-12)=>assert(Math.abs(a-b)<=tolerance,`${a} != ${b}`);
const check=(name,fn)=>{fn();checks.push({name,passed:true});};
const config={...S.defaultStudyConfig(),sampleCount:20,assessmentCount:20,intervalReplicates:50,seed:13579};
const nSchemes=config.schemes.length;
const calls=[];
const recordedPhysics={...P,solve:c=>{calls.push(clone(c));return P.solve(c);}};
const inputBefore=JSON.stringify(config);
const result=createStudy(recordedPhysics,D).runStudy(config);

check('Default candidates retain M0-M3 and add the distinct gravity-plus-Magnus combination',()=>{
 assert.deepEqual(config.schemes.map(s=>s.id),['M0','M1','M2','M3','G_M']);
 assert.deepEqual(config.schemes.map(s=>s.preset),['gravity','constantAccel','drag','magnus','fields']);
 assert.deepEqual(config.schemes.find(s=>s.id==='G_M').forceMask,{drag:false,magnus:true});
 assert.equal(S.generateSchemes({includeWind:true}).length,8);
 assert.equal(S.generateSchemes({includeWind:true,wind:[0,0,0]}).length,5);
 assert.equal(S.generateSchemes({allowDrag:false}).length,3);
 assert.equal(S.generateSchemes({allowDrag:false,allowMagnus:false,includeConstantAccel:false}).length,1);
});
check('Every pilot trial feeds exactly the same perturbed two-frame observations to every candidate',()=>{
 assert.equal(result.trials.length,20);assert.equal(calls.length,nSchemes+nSchemes*20+result.assessment.trialCount*result.assessment.models.length);
 result.trials.forEach((trial,index)=>{
  for(let k=0;k<nSchemes;k++)assert.deepEqual(calls[nSchemes+nSchemes*index+k].observations,trial.observations);
  for(let j=0;j<2;j++)for(let k=0;k<3;k++)close(trial.observations[j].position[k]-config.base.observations[j].position[k],trial.noise_m[3*j+k]);
  result.schemes.forEach(s=>assert.equal(s.samples[index].trialSeed,trial.seed));
 });
 assert.equal(JSON.stringify(config),inputBefore);
 assert(result.trials.some(t=>t.noise_m.some(x=>x!==0)));
});
check('Independent assessment shares its own trials and freezes the pilot candidate selection',()=>{
 const a=result.assessment,start=nSchemes*(1+config.sampleCount),n=a.models.length;
 assert(a.enabled);assert(a.selectionFrozen);assert.equal(a.trialCount,config.assessmentCount);
 assert.notEqual(a.seed,config.seed);assert.notDeepEqual(a.trials[0].noise_m,result.trials[0].noise_m);
 assert.deepEqual(a.models.map(s=>s.id),result.schemes.filter(s=>result.decision.candidateIds.includes(s.id)).map(s=>s.id));
 a.trials.forEach((trial,index)=>{for(let k=0;k<n;k++)assert.deepEqual(calls[start+n*index+k].observations,trial.observations);});
 a.models.forEach(s=>{assert.equal(s.samples.length,a.trialCount);assert.equal(s.density.planeArrivalMass.nTotal,a.trialCount);});
 const pilotOnly=S.runStudy({...config,assessmentEnabled:false});assert.deepEqual(pilotOnly.decision,result.decision);
 assert.equal(pilotOnly.assessment.trialCount,0);assert.deepEqual(pilotOnly.assessment.models,[]);
});
check('All attempted samples remain in state counts and observed-density denominators',()=>{
 result.schemes.forEach(s=>{
  assert.equal(s.samples.length,20);assert.equal(s.summary.total,20);
  assert.equal(Object.values(s.summary.statusCounts).reduce((a,b)=>a+b,0),20);
  assert.equal(s.summary.arrived,s.samples.filter(x=>x.status==='crossed_goal_plane').length);
  assert.equal(s.density.planeArrivalMass.nTotal,20);assert.equal(s.density.planeArrivalMass.nArrival,s.summary.arrived);
 });
});
check('Zero declared noise produces deterministic outcomes with no density or recommendation',()=>{
 const r=S.runStudy({...config,sigma:0});
 assert(r.trials.every(t=>t.noise_m.every(x=>x===0)));
 r.schemes.forEach(s=>{assert.equal(s.density.status,'degenerate');assert.equal(s.summary.eligible,false);});
 assert.equal(r.decision.preferredId,null);assert.equal(r.decision.uniqueSupported,false);
});
check('Extreme incompatible observations never produce a forced winner',()=>{
 const base=clone(config.base);base.observations[0].position[0]=20;base.observations[1].position[0]=40;
 const r=S.runStudy({...config,base});
 assert.equal(r.decision.preferredId,null);assert.equal(r.decision.uniqueSupported,false);
 r.schemes.forEach(s=>{assert.equal(s.summary.eligible,false);assert(s.summary.exclusionReasons.length);assert.equal(s.samples.length,20);});
});
let identical;
check('Variable candidate sets with identical presets and specified wind remain identical and ambiguous',()=>{
 const schemes=['A','B','C','D'].map(id=>({id,preset:'wind',label:id}));
 identical=S.runStudy({...config,schemes,wind:[0,2,0]});
 const first=identical.schemes[0];
 for(const s of identical.schemes.slice(1)){
  assert.deepEqual(s.samples,first.samples);assert.deepEqual(s.nominal.theta,first.nominal.theta);
  assert.deepEqual(s.summary,first.summary);assert.deepEqual(s.density,first.density);
 }
 assert.equal(identical.decision.uniqueSupported,false);
 if(first.summary.eligible)assert.equal(identical.decision.candidateIds.length,4);
 identical.comparison.pairs.forEach(p=>{close(p.nominalSeparation_m,0);close(p.pairedPredictionSeparationRMS_m,0);});
});
check('M1 stores and independently reconstructs all six fitted kinematic quantities for every trial',()=>{
 const m1=result.schemes.find(s=>s.id==='M1');assert(m1);
 [30.4,9.4,6.9].forEach((v,j)=>close(m1.nominal.theta[j],v));
 [-16,-16,-12].forEach((a,j)=>close(m1.nominal.fit.constantAcceleration[j],a));
 close(m1.nominal.prediction.tGoal,.405008361227,1e-10);
 const h=.05,r0=config.base.initial;
 m1.samples.forEach((sample,index)=>{
  const [r1,r2]=result.trials[index].observations.map(o=>o.position);
  const velocity=r0.map((x,j)=>(4*r1[j]-r2[j]-3*x)/(2*h));
  const acceleration=r0.map((x,j)=>(r2[j]-2*r1[j]+x)/(h*h));
  assert(sample.constantAcceleration);velocity.forEach((v,j)=>close(sample.theta[j],v,1e-9));
  acceleration.forEach((a,j)=>close(sample.constantAcceleration[j],a,1e-8));
  if(sample.status==='crossed_goal_plane'){
   const t=sample.tGoal;close(r0[0]+velocity[0]*t+.5*acceleration[0]*t*t,config.base.goalX,1e-8);
   close(r0[1]+velocity[1]*t+.5*acceleration[1]*t*t,sample.yGoal,1e-8);
   close(r0[2]+velocity[2]*t+.5*acceleration[2]*t*t,sample.zGoal,1e-8);
  }
 });
 const arrivals=m1.samples.filter(s=>s.status==='crossed_goal_plane'),n=arrivals.length;
 const meanY=arrivals.reduce((sum,s)=>sum+s.yGoal,0)/n,meanZ=arrivals.reduce((sum,s)=>sum+s.zGoal,0)/n;
 const empiricalRms=Math.sqrt(arrivals.reduce((sum,s)=>sum+(s.yGoal-meanY)**2+(s.zGoal-meanZ)**2,0)/(n-1));
 close(m1.summary.stabilityRMS_m,empiricalRms,1e-12);
 assert(m1.nominal.fit.rmse<1e-12);assert.equal(m1.nominal.fit.diagnostics.freeParameters.length,6);
});
check('Ranking and numerical dispersion do not depend on KDE bandwidth',()=>{
 const wide=S.runStudy({...config,bandwidth:{sigma_h:1.5,minBandwidth:.001}});
 assert.deepEqual(wide.decision,result.decision);
 wide.schemes.forEach((s,i)=>{
  assert.deepEqual(s.samples,result.schemes[i].samples);
  close(s.summary.stabilityRMS_m,result.schemes[i].summary.stabilityRMS_m);
  assert.deepEqual(s.summary.stabilityInterval_m,result.schemes[i].summary.stabilityInterval_m);
  if(s.density.status==='ok')assert.equal(s.density.h,1.5);
 });
});
check('Complete JSON export restores the sample ledger and every available density function',()=>{
 const restored=JSON.parse(JSON.stringify(result));
 assert.deepEqual(restored.trials,result.trials);assert.deepEqual(restored.decision,result.decision);
 restored.schemes.forEach((s,i)=>{
  assert.deepEqual(s.samples,result.schemes[i].samples);
  if(s.density.status==='ok'){
   assert.deepEqual(D.evaluateDensity(s.density,1,1),D.evaluateDensity(result.schemes[i].density,1,1));
   assert.deepEqual(D.integrateRectangle(s.density,{yMin:-3.66,yMax:3.66,zMin:.11,zMax:2.44}),D.integrateRectangle(result.schemes[i].density,{yMin:-3.66,yMax:3.66,zMin:.11,zMax:2.44}));
  }
 });
});
check('All failed fits remain in the full denominator and cannot become an empty conditional winner',()=>{
 const F={...P,solve:()=>({theta:null,fit:{ok:false,rmse:null,boundaryHit:false},prediction:{status:'fit_failed'},issues:['Synthetic optimizer failure for ledger test']})};
 const r=createStudy(F,D).runStudy(config);
 r.schemes.forEach(s=>{
  assert.equal(s.summary.total,20);assert.equal(s.summary.arrived,0);assert.equal(s.summary.statusCounts.fit_failed,20);
  assert.equal(s.summary.numericalFailureFraction,1);assert.equal(s.density.status,'degenerate');
  assert.equal(s.density.planeArrivalMass.nTotal,20);assert.equal(s.density.planeArrivalMass.fraction,0);
 });
 assert.equal(r.decision.preferredId,null);assert.equal(r.decision.uniqueSupported,false);
});
check('Thrown solver failures and each nonarrival status are retained separately',()=>{
 let ncall=0;
 const states=['crossed_goal_plane','ground_before_goal','no_crossing_within_horizon','invalid_input','fit_failed','integration_failed','outside_model_domain'];
 const F={...P,solve:c=>{
  const call=ncall++;if(call<nSchemes)return P.solve(c);
  const i=Math.floor((call-nSchemes)/nSchemes),status=states[i%states.length];
  if(status==='integration_failed')throw Error('Synthetic integration failure');
  return {theta:[1,2,3],fit:{ok:status==='crossed_goal_plane',rmse:0,boundaryHit:i%2===0},prediction:{status,yGoal:status==='crossed_goal_plane'?i/10:null,zGoal:status==='crossed_goal_plane'?1:null,tGoal:status==='crossed_goal_plane'?.4:null},issues:[]};
 }};
 const r=createStudy(F,D).runStudy(config);
 r.schemes.forEach(s=>{
  assert.equal(s.samples.length,20);assert.equal(Object.values(s.summary.statusCounts).reduce((a,b)=>a+b,0),20);
  states.forEach(state=>assert(s.summary.statusCounts[state]>0));
  assert.equal(s.summary.arrived,3);assert.equal(s.density.planeArrivalMass.nTotal,20);
  assert.equal(s.summary.eligible,false);
 });
});
check('Invalid smoothing floor is rejected before producing a recommendation',()=>{
 assert(S.validate({...config,bandwidth:{sigma_h:null,minBandwidth:0}}).length>0);
});
check('Configured minimum sample threshold prevents a small-sample unique recommendation',()=>{
 const r=S.runStudy({...config,sampleCount:2});
 assert.equal(r.decision.uniqueSupported,false);assert.equal(r.decision.preferredId,null);
 r.schemes.forEach(s=>assert(s.summary.exclusionReasons.includes('可用样本不足')));
});
check('All three factories serialize together and preserve the paired random observations',()=>{
 const factory='(()=>{const P=('+createPhysics.toString()+')();const D=('+createDensity.toString()+')();return ('+createStudy.toString()+')(P,D);})()';
 const worker=vm.runInNewContext(factory);
 assert.deepEqual(JSON.parse(JSON.stringify(worker.makeTrials(clone(config)))),result.trials);
});
const report={passed:true,checks:checks.length,node:process.version,studyVersion:S.VERSION,densityVersion:D.version,tests:checks,
 sampleCount:config.sampleCount,seed:config.seed,defaultDecision:result.decision,
 scope:'Checks software comparison and all-trial bookkeeping, not physical validity or a uniquely correct force model.'};
if(process.argv.includes('--write-report'))fs.writeFileSync(path.join(__dirname,'study.validation.json'),JSON.stringify(report,null,2));
console.log(JSON.stringify({passed:true,checks:checks.length,defaultDecision:result.decision.status},null,2));
