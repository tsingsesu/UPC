/* Run: node src/web/physics.test.cjs [--write-report]
   No writes unless explicitly requested. Tests use independent MATLAB fixtures,
   exact solutions and accounting invariants; no browser or network required. */
'use strict';
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const P=require('./physics.js');
const checks=[],began=Date.now();
const maxDiff=(a,b)=>Math.max(...a.map((x,i)=>Math.abs(x-b[i])));
function near(a,b,tol,label){assert.ok(Number.isFinite(a)&&Math.abs(a-b)<=tol,`${label}: ${a} vs ${b}, tol ${tol}`);}
function test(name,fn){const t=Date.now();try{const evidence=fn();checks.push({name,passed:true,elapsedMs:Date.now()-t,evidence:evidence||null});}catch(e){checks.push({name,passed:false,error:e.stack});}}
const fixture=JSON.parse(fs.readFileSync(path.join(__dirname,'../../data/验证基准/浏览器对照.json'),'utf8'));
let nominal;
test('V08 MATLAB eight-case cross implementation',()=>{
  const comparisons=[];for(const f of fixture.cases){const r=P.solve(f.browserConfig);
    assert.equal(r.prediction.status,fixture.statusMapping[f.prediction.status]||f.prediction.status,f.id);
    let eventError=null;if(r.prediction.status==='crossed_goal_plane'){
      eventError={y:Math.abs(r.prediction.yGoal-f.prediction.yGoal),z:Math.abs(r.prediction.zGoal-f.prediction.zGoal),t:Math.abs(r.prediction.tGoal-f.prediction.tGoal)};
      near(r.prediction.yGoal,f.prediction.yGoal,1e-4,f.id+' y');near(r.prediction.zGoal,f.prediction.zGoal,1e-4,f.id+' z');near(r.prediction.tGoal,f.prediction.tGoal,1e-5,f.id+' t');}
    const positions=P.predictPositions(r.config,r.theta,f.sampleTimes);
    const positionError=Math.max(...positions.map((p,i)=>maxDiff(p,f.samplePositions[i])));assert.ok(positionError<1e-4,f.id+' positions');
    const observationError=Math.max(...r.predictedObservations.map((p,i)=>maxDiff(p,f.observationsPredicted[i])));assert.ok(observationError<1e-5,f.id+' observations');
    comparisons.push({id:f.id,eventError,positionError,observationError});if(f.id==='inverse_M3')nominal=r;
  }return {matlabVersion:fixture.matlabVersion,comparisons};
});
test('V08 four force cases, SI and sign',()=>{
  return fixture.forceCases.map(f=>{const a=P.evaluateForces(f.t,f.state,f.browserConfig,f.theta);let error=0;
    for(const k of ['acceleration','gravity','drag','magnus','total'])error=Math.max(error,maxDiff(a[k],f.expected[k==='total'?'totalForce':k]));
    assert.ok(error<1e-8);return {id:f.id,maxError:error};});
});
test('V01 gravity exact trajectory and event',()=>{
  const c={...P.defaultConfig(),mode:'forward',model:'M0',velocity:[25,4,6],h:.003},r=P.solve(c),t=11/25;
  near(r.prediction.tGoal,t,1e-9,'gravity time');near(r.prediction.yGoal,4*t,1e-9,'gravity y');near(r.prediction.zGoal,.11+6*t-c.g*t*t/2,1e-9,'gravity z');
  let error=0;for(const q of r.trajectory)error=Math.max(error,maxDiff(q.position,[25*q.t,4*q.t,.11+6*q.t-c.g*q.t*q.t/2]));assert.ok(error<1e-9);return {maxPositionError:error};
});
test('V02 one-dimensional quadratic drag exact solution',()=>{
  const c={...P.defaultConfig(),mode:'forward',model:'M2',g:0,initial:[0,0,1],velocity:[18,0,0],Cd:.25};const r=P.solve(c),k=c.rho*Math.PI*c.R*c.R*c.Cd/(2*c.m);
  let error=0;for(const q of r.trajectory){error=Math.max(error,Math.abs(q.position[0]-Math.log(1+k*18*q.t)/k),Math.abs(q.velocity[0]-18/(1+k*18*q.t)));}assert.ok(error<1e-8);return {maxError:error};
});
test('V03 energy budget matches integrated aerodynamic power',()=>{
  const c={...P.defaultConfig(),mode:'forward',velocity:[31,8,7],Cd:.35,b:[0,-.1,-.3]},r=P.solve(c);let work=0;
  const energy=q=>.5*c.m*q.velocity.reduce((s,x)=>s+x*x,0)+c.m*c.g*q.position[2];
  const power=q=>q.forces.drag.reduce((s,x,i)=>s+x*q.velocity[i],0)+q.forces.magnus.reduce((s,x,i)=>s+x*q.velocity[i],0);
  for(let i=1;i<r.trajectory.length;i++){const a=r.trajectory[i-1],b=r.trajectory[i];work+=(power(a)+power(b))*.5*(b.t-a.t);assert.ok(energy(b)<=energy(a)+1e-8);}
  const error=Math.abs(energy(r.trajectory.at(-1))-energy(r.trajectory[0])-work);assert.ok(error<1e-4);return {energyWorkError_J:error};
});
test('V04/V05 Magnus perpendicular, zero relative velocity and parallel b',()=>{
  const c={...P.defaultConfig(),mode:'forward',b:[.2,-.3,.4]},th=[31,8,7,.3,-.3,.4],f=P.evaluateForces(0,[0,0,.11,31,8,7],c,th);
  near(f.magnus.reduce((s,x,i)=>s+x*f.relativeAirVelocity[i],0),0,1e-12,'Magnus power');
  const zero=P.evaluateForces(0,[0,0,1,0,0,0],c,th);assert.ok(zero.drag.every(x=>x===0));assert.ok(zero.magnus.every(x=>x===0));
  assert.equal(zero.CL_effective,null);const par=P.evaluateForces(0,[0,0,1,1,2,3],{...c,b:[1,2,3]},[1,2,3,.3,2,3]);assert.deepEqual(par.magnus,[0,0,0]);
});
test('V06 event states and initial takeoff',()=>{
  const d={...P.defaultConfig(),mode:'forward',model:'M0'};
  const ground=P.solve({...d,velocity:[1,0,-1]});assert.equal(ground.prediction.status,'ground_before_goal');assert.equal(ground.prediction.eventTime,0);assert.equal(ground.prediction.tGoal,null);assert.equal(ground.predictedObservations[0],null);
  assert.equal(P.solve({...d,velocity:[30,0,8]}).prediction.status,'crossed_goal_plane');
  const falling=P.solve({...d,velocity:[1,0,1]});assert.equal(falling.prediction.status,'ground_before_goal');near(falling.prediction.eventTime,2/d.g,1e-9,'return ground');
  assert.equal(P.solve({...d,initial:[0,0,5],velocity:[0,0,0],g:0}).prediction.status,'no_crossing_within_horizon');
  assert.equal(P.solve({...d,R:.12}).prediction.status,'invalid_input');
  assert.equal(P.integrateTrajectory(d,[1000,0,1,.3,0,0]).prediction.status,'outside_model_domain');
  assert.equal(P.integrateTrajectory(d,[1,2]).prediction.status,'invalid_input');
  assert.equal(P.solve({...d,h:NaN}).prediction.status,'invalid_input');
  assert.equal(P.solve({...d,unknownTypo:1}).prediction.status,'invalid_input');
  const advanced={...d,omega:[0,0,50],quaternion:[1,0,0,0],spinDecayTau:2,kappa:2/3,dragModel:'constant',ballType:'generic',liftModel:'linear',liftSlope:1,mu:1.8e-5,windShear:[0,0,0],gustAmplitude:[0,0,0],gustPeriod:1};
  const r=P.solve(advanced);assert.ok(r.theta);assert.ok(r.parameterRoles.inactive.includes('omega'));assert.deepEqual(r.config.omega,[0,0,50]);
  assert.equal(P.solve({...d,metadata:{note:'preserved and not applied'}}).config.metadata.note,'preserved and not applied');
  assert.throws(()=>P.predictPositions(d,[1,2,3,.3,0,0],[Infinity]));
});
test('V07 three independent integration step sizes',()=>{
  nominal=nominal||P.solve(P.defaultConfig());const c={...P.defaultConfig(),mode:'forward',velocity:nominal.theta.slice(0,3),Cd:nominal.theta[3],b:[0,...nominal.theta.slice(4)]};
  const r=[.004,.002,.001].map(h=>P.solve({...c,h}).prediction),changes=r.slice(0,-1).map((a,i)=>({step:[.004,.002][i],position:Math.hypot(a.yGoal-r[i+1].yGoal,a.zGoal-r[i+1].zGoal),time:Math.abs(a.tGoal-r[i+1].tGoal)}));
  assert.ok(changes.every(e=>e.position<1e-4&&e.time<1e-5));return {changes};
});
test('M1 algebra regression and all six observation edits truly refit',()=>{
  const c=P.defaultConfig(),m1=P.solve({...c,model:'M1'});near(m1.prediction.tGoal,.4050083612,1e-9,'M1 t');near(m1.prediction.yGoal,2.4948244142,1e-9,'M1 y');near(m1.prediction.zGoal,1.9203670565,1e-9,'M1 z');
  const base=P.solve(c),changes=[];for(let j=0;j<6;j++){const x=P.defaultConfig();x.observations[Math.floor(j/3)].position[j%3]+=.005;const r=P.solve(x);assert.ok(r.fit.ok);const change=Math.hypot(r.prediction.yGoal-base.prediction.yGoal,r.prediction.zGoal-base.prediction.zGoal);assert.ok(change>1e-4);changes.push(change);}return {goalChanges_m:changes};
});
test('identifiability diagnostics and synthetic same-model recovery',()=>{
  const cases=[];for(const th of [[28,4,6,.2,.1,-.2],[35,6,8,.45,-.2,-.4],[24,3,5,.3,.2,.1]]){
    const c={...P.defaultConfig(),mode:'forward',velocity:th.slice(0,3),Cd:th[3],b:[0,...th.slice(4)]},truth=P.solve(c);
    const inverse={...c,mode:'inverse',observations:truth.predictedObservations.map((position,i)=>({t:[.05,.1][i],position}))},r=P.solve(inverse);
    assert.ok(r.fit.ok);assert.equal(r.fit.diagnostics.rank,6);assert.equal(r.fit.starts.length,4);assert.ok(r.fit.starts.every(s=>s.theta.length===6));
    const error=Math.hypot(r.prediction.yGoal-truth.prediction.yGoal,r.prediction.zGoal-truth.prediction.zGoal);assert.ok(error<1e-5);cases.push({theta:th,goalError:error});}return {cases};
});
test('V10 seeded bootstrap reproducibility and complete accounting',()=>{
  const c=P.defaultConfig(),a=P.uncertainty(c,24,17),b=P.uncertainty(c,24,17);
  assert.deepEqual(a.samples,b.samples);const sum=Object.entries(a.summary).filter(([k])=>['crossed_goal_plane','ground_before_goal','no_crossing_within_horizon','invalid_input','fit_failed','integration_failed','outside_model_domain','cancelled'].includes(k)).reduce((s,[,n])=>s+n,0);assert.equal(sum,24);
  assert.equal(a.summary.denominators.arrivalFraction,24);assert.equal(a.summary.denominators.centerInsideFractionOfArrivals,a.summary.reached);
  assert.throws(()=>P.uncertainty({...c,sigma:0},10,1));
  const fake=[{status:'fit_failed'},{status:'integration_failed'},{status:'ground_before_goal'},{status:'crossed_goal_plane',yGoal:0,zGoal:1,tGoal:.4,centerInside:true,clearanceInside:true}];const s=P.describeSamples(fake,'bootstrap');assert.equal(s.total,4);assert.equal(s.arrivalFraction,.25);assert.equal(s.centerInsideFractionOfArrivals,1);
  return {sampleCount:24,summary:a.summary};
});
test('deterministic scenario has 12 entries and no probability ellipse',()=>{
  const r=P.uncertainty(P.defaultConfig(),20,123,null,'scenario');assert.equal(r.samples.length,12);assert.equal(r.summary.conditionalEllipse,null);assert.equal(r.summary.yQuantiles,null);return {summary:r.summary};
});
test('worker serialization and JSON roundtrip',()=>{
  const factory=Function('return ('+P.createPhysics.toString()+')')(),c=JSON.parse(JSON.stringify(P.defaultConfig()));const a=P.solve(c),b=factory().solve(c);assert.deepEqual(a.theta,b.theta);assert.deepEqual(a.prediction,b.prediction);
});
test('incompatible observations distinguish optimizer success from explanation',()=>{
  const c=P.defaultConfig();c.observations[0].position[0]=20;c.observations[1].position[0]=40;
  const r=P.solve(c);assert.ok(r.fit.ok);assert.equal(r.fit.diagnostics.compatibility.flag,'incompatible_with_declared_noise');assert.ok(r.issues.length>0);assert.ok(r.fit.boundaryHit);
  const noNoise=P.solve({...c,sigma:0});assert.equal(noNoise.fit.diagnostics.compatibility.flag,'large_residual');
  return {rmse_m:r.fit.rmse,compatibility:r.fit.diagnostics.compatibility,boundaryHit:r.fit.boundaryHit};
});
const report={schemaVersion:'1.0',codeVersion:P.VERSION,scope:'Numerical verification of selected equations, not experimental confirmation.',node:process.version,platform:process.platform,checks,passed:checks.every(t=>t.passed),elapsedMs:Date.now()-began};
if(process.argv.includes('--write-report')){const dest=path.join(__dirname,'../../data/验证基准/浏览器内核测试.json');fs.writeFileSync(dest,JSON.stringify(report,null,2),'utf8');}
console.log(JSON.stringify(report,null,2));
if(!report.passed)process.exitCode=1;
