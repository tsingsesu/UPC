'use strict';
const assert=require('node:assert/strict');
const {createPhysics}=require('./physics.js');
const {createAdvanced}=require('./advanced.js');
const fs=require('node:fs'),path=require('node:path'),dir=path.join(__dirname,'..','..','data','验证基准');
const A=createAdvanced(createPhysics());
const tests=[];
const crossImplementation=[];
const n=v=>Math.hypot(...v),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const near=(a,b,tol=1e-8)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);
function run(name,f){f();tests.push(name);}
function cfg(extra={}){return {...A.defaultConfig(),model:'M4',mode:'forward',...extra};}
function state(c){return [...c.initial,...c.velocity,...c.omega,...c.quaternion];}
run('M0-M3 delegation unchanged',()=>{
  const P=createPhysics(),c=P.defaultConfig();
  const a=A.solve(c),b=P.solve(c);
  near(a.prediction.yGoal,b.prediction.yGoal,1e-12);
  near(a.prediction.zGoal,b.prediction.zGoal,1e-12);
  assert.equal(P.solve(A.defaultConfig()).prediction.status,'crossed_goal_plane');
});
run('P0 comparisons and bootstrap remain available with advanced fields',()=>{
  const c=A.defaultConfig(),comp=A.compare(c,'refit'),sample=A.uncertainty(c,3,20261107);
  assert.equal(comp.runs.length,4);assert.ok(comp.runs.every(r=>r.prediction.status==='crossed_goal_plane'));
  assert.equal(sample.samples.length,3);assert.equal(sample.summary.total,3);
});
run('M4 forward only and malformed input',()=>{
  assert.equal(A.solve(cfg({mode:'inverse'})).prediction.status,'invalid_input');
  assert.equal(A.solve(cfg({quaternion:[2,0,0,0]})).prediction.status,'invalid_input');
  assert.equal(A.solve(cfg({velocity:1})).prediction.status,'invalid_input');
  assert.throws(()=>A.uncertainty(cfg(),10),/概率推断/);
});
run('Gravity analytic with no aerodynamic terms',()=>{
  const c=cfg({velocity:[25,4,6],Cd:0,liftSlope:0});
  const r=A.solve(c),t=11/25;
  assert.equal(r.prediction.status,'crossed_goal_plane');
  near(r.prediction.tGoal,t,1e-9);
  near(r.prediction.yGoal,4*t,1e-9);
  near(r.prediction.zGoal,.11+6*t-.5*c.g*t*t,1e-9);
});
run('Constant world spin gives analytic scalar-first quaternion',()=>{
  const c=cfg({g:0,initial:[0,0,1],velocity:[10,0,0],Cd:0,liftSlope:0,omega:[0,0,50]});
  const r=A.solve(c),last=r.trajectory.at(-1),phase=50*last.t/2;
  near(n(last.q),1,1e-12);near(last.q[0],Math.cos(phase),2e-7);near(last.q[3],Math.sin(phase),2e-7);
  near(n(last.omega),50,1e-12);
});
run('Exponential spin decay, angular momentum and quaternion phase',()=>{
  const c=cfg({g:0,initial:[0,0,1],velocity:[10,0,0],Cd:0,liftSlope:0,omega:[30,0,0],spinDecayTau:.3});
  const r=A.solve(c),last=r.trajectory.at(-1),expectedOmega=30*Math.exp(-last.t/.3);
  const angle=30*.3*(1-Math.exp(-last.t/.3));
  near(last.omega[0],expectedOmega,1e-8);near(last.q[0],Math.cos(angle/2),1e-7);near(last.q[1],Math.sin(angle/2),1e-7);
  const I=c.kappa*c.m*c.R*c.R;
  near(last.forces.torque[0],-I*expectedOmega/.3,1e-10);
});
run('Magnus work is zero in air frame; zero spin and parallel axis limits',()=>{
  let c=cfg(),f=A.evaluateForces(0,state(c),c);
  near(dot(f.magnus,f.relativeAirVelocity),0,1e-10);
  c=cfg({omega:[0,0,0]});f=A.evaluateForces(0,state(c),c);near(n(f.magnus),0,1e-12);
  c=cfg({velocity:[30,0,0],omega:[50,0,0]});f=A.evaluateForces(0,state(c),c);near(n(f.magnus),0,1e-12);
});
run('Zero relative speed produces no aerodynamic force or NaN',()=>{
  const c=cfg({velocity:[20,1,2],wind:[20,1,2]});const f=A.evaluateForces(0,state(c),c);
  near(n(f.drag),0,1e-12);near(n(f.magnus),0,1e-12);
  assert.equal(f.S,null);assert.ok(f.acceleration.every(Number.isFinite));
});
run('Preset shear and sinusoidal gust feed relative velocity',()=>{
  const c=cfg({initial:[0,0,1],velocity:[30,0,0],wind:[1,2,3],windShear:[0,2,0],gustAmplitude:[0,0,4],gustPeriod:2});
  const f=A.evaluateForces(.5,state(c),c);
  assert.deepEqual(f.airVelocity,[1,4,7]);assert.deepEqual(f.relativeAirVelocity,[29,-4,-7]);
});
run('Paper2018 Table 1 equation 5 independent numeric check',()=>{
  for(const ballType of ['Tango12','Teamgeist','Jabulani','Brazuca']){
    const c=cfg({dragModel:'paper2018',ballType,velocity:[30,0,0],omega:[0,0,0]});
    const [a,vc,vs,bmin,bmax,vmin,vmax]=A.advancedDragTable()[ballType];
    const expected=(a-bmin)/(1+Math.exp((30-vc)/vs))+bmin+(30-vmin)/(1+Math.exp(-(30-vmin)/vs))*(bmax-bmin)/(vmax-vmin);
    const f=A.evaluateForces(0,state(c),c);near(f.Cd,expected,1e-12);assert.equal(f.domainStatus,'ok');
  }
});
run('Teamgeist spin correction and lift equation 9',()=>{
  const c=cfg({dragModel:'paper2018',liftModel:'teamgeist2018',ballType:'Teamgeist',velocity:[30,0,0],omega:[0,0,50]});
  const f=A.evaluateForces(0,state(c),c),S=.11*50/30;
  const cd0=(.4927-.144)/(1+Math.exp((30-12.58)/1.071))+.144+(30-23.17)/(1+Math.exp(-(30-23.17)/1.071))*.01/(35-23.17);
  near(f.Cd,cd0+.514*S,1e-12);
  near(f.CL,1.15*S**.83*(.4927-Math.min(f.Cd,.4927))/(.4927-(.144+.514*.19)),1e-12);
  near(dot(f.magnus,f.relativeAirVelocity),0,1e-10);
});
run('Domain exits stop with explicit status and no claimed goal',()=>{
  const high=A.solve(cfg({dragModel:'paper2018',velocity:[36,0,1]}));
  assert.equal(high.prediction.status,'outside_model_domain');assert.equal(high.prediction.eventTime,0);assert.equal(high.prediction.tGoal,null);
  const spin=A.solve(cfg({dragModel:'paper2018',liftModel:'teamgeist2018',omega:[0,0,150]}));
  assert.equal(spin.prediction.status,'outside_model_domain');
  const slowing=A.solve(cfg({initial:[0,0,1],velocity:[7.1,0,0],g:0,liftSlope:0,dragModel:'paper2018',goalX:100}));
  assert.equal(slowing.prediction.status,'outside_model_domain');assert.ok(slowing.prediction.eventTime>0);
  near(slowing.trajectory.at(-1).forces.airSpeed,7,1e-8);
});
run('First-flight ground event and unavailable post-event observations',()=>{
  const r=A.solve(cfg({velocity:[25,0,.1],Cd:0,liftSlope:0}));
  assert.equal(r.prediction.status,'ground_before_goal');near(r.prediction.eventTime,.2/9.80665,1e-8);
  assert.equal(r.prediction.tGoal,null);assert.deepEqual(r.predictedObservations,[null,null]);
});
run('Mechanical plus rotational energy closes drag and torque work',()=>{
  const c=cfg({initial:[0,0,1],velocity:[25,0,3],omega:[0,0,50],spinDecayTau:.5,h:.0005});
  const r=A.solve(c),rows=r.trajectory,I=c.kappa*c.m*c.R*c.R;
  const energy=x=>.5*c.m*dot(x.velocity,x.velocity)+c.m*c.g*x.position[2]+.5*I*dot(x.omega,x.omega);
  const power=x=>dot(x.forces.drag,x.velocity)+dot(x.forces.magnus,x.velocity)+dot(x.torque,x.omega);
  let work=0;for(let i=1;i<rows.length;i++)work+=(rows[i].t-rows[i-1].t)*(power(rows[i])+power(rows[i-1]))/2;
  near(energy(rows.at(-1))-energy(rows[0]),work,1e-4);
});
run('Worker serialization retains complete advanced implementation',()=>{
  const factory=eval('('+createAdvanced.toString()+')'),B=factory(createPhysics());
  const c=cfg();near(B.solve(c).prediction.yGoal,A.solve(c).prediction.yGoal,1e-12);
});
run('Independent MATLAB ode45 13-state comparison across five scenarios',()=>{
  const reference=JSON.parse(fs.readFileSync(path.join(dir,'M4跨实现对照.json'),'utf8'));
  for(const test of reference.cases){
    const result=A.solve(test.config),expected=test.prediction;
    assert.equal(result.prediction.status,expected.status,test.id);
    const timeError=Math.abs(result.prediction.eventTime-expected.eventTime);assert.ok(timeError<1e-5,test.id+' event time');
    let positionError=0,velocityError=0,omegaError=0,quaternionError=0;
    const sampleTimes=Array.isArray(test.sampleTimes)?test.sampleTimes:[test.sampleTimes];
    const states=Array.isArray(test.states[0])?test.states:[test.states];
    function compareRow(actual,truth){
      positionError=Math.max(positionError,...actual.position.map((x,i)=>Math.abs(x-truth[i])));
      velocityError=Math.max(velocityError,...actual.velocity.map((x,i)=>Math.abs(x-truth[3+i])));
      omegaError=Math.max(omegaError,...actual.omega.map((x,i)=>Math.abs(x-truth[6+i])));
      const qm=Math.max(...actual.q.map((x,i)=>Math.abs(x-truth[9+i]))),qp=Math.max(...actual.q.map((x,i)=>Math.abs(x+truth[9+i])));
      quaternionError=Math.max(quaternionError,Math.min(qm,qp));
    }
    for(let i=0;i<sampleTimes.length;i++){
      const target=sampleTimes[i],row=result.trajectory.reduce((best,r)=>Math.abs(r.t-target)<Math.abs(best.t-target)?r:best);
      assert.ok(Math.abs(row.t-target)<1e-10,test.id+' sample time');compareRow(row,states[i]);
    }
    compareRow(result.trajectory.at(-1),test.finalState);
    assert.ok(positionError<1e-4,test.id+' position '+positionError);
    assert.ok(velocityError<1e-4,test.id+' velocity '+velocityError);
    assert.ok(omegaError<1e-5,test.id+' omega '+omegaError);
    assert.ok(quaternionError<1e-5,test.id+' quaternion '+quaternionError);
    crossImplementation.push({id:test.id,status:expected.status,timeError_s:timeError,positionError_m:positionError,
      velocityError_m_s:velocityError,omegaError_rad_s:omegaError,quaternionError,referenceMethod:reference.method});
  }
});
const report={passed:true,count:tests.length,tests,crossImplementation,nodeVersion:process.version,physicalValidation:'not independently experimentally validated'};
fs.mkdirSync(dir,{recursive:true});fs.writeFileSync(path.join(dir,'高级模型测试.json'),JSON.stringify(report,null,2),'utf8');
console.log(JSON.stringify(report,null,2));
