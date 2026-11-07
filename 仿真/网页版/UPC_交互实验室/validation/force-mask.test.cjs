'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const source=path.resolve(__dirname,'../physics.js');
const {createPhysics}=require(source);
const P=createPhysics(),legacy=JSON.parse(fs.readFileSync(path.join(__dirname,'legacy_nominal_reference.json'),'utf8')),checks=[];
const clone=x=>JSON.parse(JSON.stringify(x));
const norm=x=>Math.hypot(...x),dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
const maxdiff=(a,b)=>Math.max(...a.map((x,i)=>Math.abs(x-b[i])));
const close=(a,b,t=1e-10)=>assert(Math.abs(a-b)<t,`${a} != ${b}`);
const check=(name,f)=>{f();checks.push({name,passed:true});};
const regression=[];
check('Original M0-M3 official nominal trajectories remain numerically unchanged',()=>{
 for(const model of ['M0','M1','M2','M3']){
  const a=legacy.cases[model],b=P.solve({...P.defaultConfig(),model});
  assert.equal(b.fit.ok,true);assert.equal(b.prediction.status,a.prediction.status);
  const error=maxdiff([b.prediction.yGoal,b.prediction.zGoal,b.prediction.tGoal],[a.prediction.yGoal,a.prediction.zGoal,a.prediction.tGoal]);
  assert(error<1e-10);assert(maxdiff(a.theta,b.theta)<1e-9);regression.push({model,maximumGoalDifference:error});
 }
});
check('Legacy JSON without forceMask defaults to both physical aerodynamic terms enabled',()=>{
 const r=P.solve(legacy.defaultConfig),baseline=P.solve(P.defaultConfig());
 assert.deepEqual(r.theta,baseline.theta);assert.deepEqual(r.config.forceMask,{drag:true,magnus:true});
});
check('Mask changes active fit dimensions and locks disabled coefficients at zero',()=>{
 for(const [drag,magnus,count] of [[true,true,6],[true,false,4],[false,true,5],[false,false,3]]){
  const r=P.solve({...P.defaultConfig(),forceMask:{drag,magnus}});
  assert(r.fit.ok);assert.equal(r.fit.diagnostics.freeParameters.length,count);
  assert.equal(r.fit.diagnostics.jacobian[0].length,count);
  if(!drag){assert.equal(r.theta[3],0);assert.equal(r.config.Cd,0);assert(!r.fit.diagnostics.freeParameters.includes('Cd'));}
  if(!magnus){assert.deepEqual(r.theta.slice(4),[0,0]);assert.deepEqual(r.config.b,[0,0,0]);assert(!r.fit.diagnostics.freeParameters.includes('by'));}
 }
});
check('Drag-only M3 matches M2 inverse and both-off M3 matches the M0 optimum',()=>{
 for(const [mask,model] of [[{drag:true,magnus:false},'M2'],[{drag:false,magnus:false},'M0']]){
  const a=P.solve({...P.defaultConfig(),forceMask:mask}),b=P.solve({...P.defaultConfig(),model});
  assert(maxdiff(a.theta,b.theta)<1e-8);close(a.fit.rmse,b.fit.rmse,1e-12);
  close(a.prediction.yGoal,b.prediction.yGoal,1e-8);close(a.prediction.zGoal,b.prediction.zGoal,1e-8);close(a.prediction.tGoal,b.prediction.tGoal,1e-9);
 }
});
check('Forward masks disable actual forces, including bx, and preserve input objects',()=>{
 for(const mask of [{drag:false,magnus:true},{drag:true,magnus:false},{drag:false,magnus:false}]){
  const c={...P.defaultConfig(),mode:'forward',velocity:[30,5,6],Cd:.6,b:[.2,.15,-.25],forceMask:mask},before=JSON.stringify(c),r=P.solve(c);
  assert.equal(JSON.stringify(c),before);assert.equal(r.prediction.status,'crossed_goal_plane');
  for(const row of r.trajectory){
   if(!mask.drag){close(norm(row.forces.drag),0);assert.equal(row.forces.Cd,0);}
   if(!mask.magnus){close(norm(row.forces.magnus),0);assert.deepEqual(row.forces.effectiveB,[0,0,0]);}
   const expected=row.forces.gravity.map((x,j)=>x+row.forces.drag[j]+row.forces.magnus[j]);
   assert(maxdiff(expected,row.forces.total)<1e-12);
  }
 }
});
function analyticMagnus(c,velocity,b,t){
 const B=norm(b),axis=b.map(x=>x/B),w=velocity.map((x,i)=>x-c.wind[i]),g=[0,0,-c.g];
 const wp=axis.map(x=>x*dot(axis,w)),wg=w.map((x,i)=>x-wp[i]),gp=axis.map(x=>x*dot(axis,g)),gg=g.map((x,i)=>x-gp[i]);
 const nw=cross(axis,w),ng=cross(axis,g),sin=Math.sin(B*t),cos=Math.cos(B*t);
 return c.initial.map((x,i)=>x+c.wind[i]*t+wp[i]*t+sin/B*wg[i]+(1-cos)/B*nw[i]+.5*gp[i]*t*t+(1-cos)/(B*B)*gg[i]+(t/B-sin/(B*B))*ng[i]);
}
check('Magnus-only trajectories agree with independent constant-b matrix-exponential formula',()=>{
 const c={...P.defaultConfig(),mode:'forward',forceMask:{drag:false,magnus:true},wind:[2,-3,.5],b:[.2,.15,-.25]},th=[30,5,6,1.5,.15,-.25],times=[.05,.1,.3];
 const predicted=P.predictPositions(c,th,times);
 times.forEach((t,i)=>assert(maxdiff(predicted[i],analyticMagnus(c,th.slice(0,3),c.b,t))<2e-10));
 const f=P.evaluateForces(.1,[...predicted[1],...th.slice(0,3)],c,th);
 close(dot(f.magnus,f.relativeAirVelocity),0,1e-12);assert.equal(f.Cd,0);
});
let syntheticError;
check('Magnus-only five-parameter inverse recovers independent analytic observations with specified wind',()=>{
 const c={...P.defaultConfig(),forceMask:{drag:false,magnus:true},wind:[2,-3,.5]},truth=[30,5,6,0,.15,-.25],b=[0,truth[4],truth[5]];
 c.observations=[.05,.1].map(t=>({t,position:analyticMagnus(c,truth.slice(0,3),b,t)}));
 const r=P.solve(c);assert(r.fit.ok);syntheticError=maxdiff(r.theta,truth);assert(syntheticError<2e-8);assert(r.fit.rmse<1e-11);
 assert.deepEqual(r.fit.diagnostics.freeParameters,['vx','vy','vz','by','bz']);assert.equal(r.theta[3],0);
});
check('Unused aerodynamic input values do not leak into inverse results',()=>{
 const a=P.solve({...P.defaultConfig(),forceMask:{drag:false,magnus:true},Cd:.05});
 const b=P.solve({...P.defaultConfig(),forceMask:{drag:false,magnus:true},Cd:1.5});
 assert.deepEqual(a.theta,b.theta);assert.deepEqual(a.prediction,b.prediction);
 const c=P.solve({...P.defaultConfig(),forceMask:{drag:true,magnus:false},b:[1,.8,-.7]});
 const d=P.solve({...P.defaultConfig(),forceMask:{drag:true,magnus:false},b:[0,0,0]});
 assert.deepEqual(c.theta,d.theta);assert.deepEqual(c.prediction,d.prediction);
});
check('M0-M2 retain their model equations and explicitly mark mask inactive',()=>{
 for(const model of ['M0','M1','M2']){
  const a=P.solve({...P.defaultConfig(),model,forceMask:{drag:false,magnus:false}}),b=P.solve({...P.defaultConfig(),model});
  assert.deepEqual(a.theta,b.theta);assert.deepEqual(a.prediction,b.prediction);assert(a.parameterRoles.inactive.includes('forceMask'));
 }
});
check('Invalid mask or explicit nonfinite theta fails rather than switching terms silently',()=>{
 for(const forceMask of [null,{drag:0,magnus:true},{drag:true},{drag:true,magnus:true,gravity:false}])assert.equal(P.solve({...P.defaultConfig(),forceMask}).prediction.status,'invalid_input');
 const c={...P.defaultConfig(),forceMask:{drag:false,magnus:false}};
 assert.equal(P.integrateTrajectory(c,[]).prediction.status,'invalid_input');
 assert.equal(P.integrateTrajectory(c,[30,5,6,NaN,0,0]).prediction.status,'invalid_input');
});
check('M1 interpolation remains correct and rejection diagnostics distinguish domain from contact',()=>{
 const c={...P.defaultConfig(),model:'M1'},r=P.solve(c);
 assert(maxdiff(r.theta.slice(0,3),[30.4,9.4,6.9])<1e-12);assert(maxdiff(r.fit.constantAcceleration,[-16,-16,-12])<1e-11);
 const v=[30,0,6],a=[-600,0,-10];
 c.observations=[.05,.1].map(t=>({t,position:c.initial.map((x,j)=>x+v[j]*t+.5*a[j]*t*t)}));
 const outside=P.solve(c);assert.equal(outside.fit.ok,false);assert.equal(outside.fit.diagnostics.rejectionReason,'invalid_input');
 assert(outside.fit.diagnostics.issues.some(x=>x.includes('500')));assert(!outside.fit.diagnostics.message.includes('观测前触地'));
 const speed={...P.defaultConfig(),model:'M1'};speed.observations=[.05,.1].map(t=>({t,position:[110*t,0,.11+6*t-5*t*t]}));
 const fast=P.solve(speed);assert.equal(fast.fit.ok,false);assert.equal(fast.fit.diagnostics.rejectionReason,'outside_model_domain');
});
const report={passed:true,checks:checks.length,version:P.VERSION,node:process.version,tests:checks,regression,
 independentMagnusSyntheticMaximumParameterError:syntheticError,source:'../physics.js',sourceSHA256:crypto.createHash('sha256').update(fs.readFileSync(source)).digest('hex'),
 legacyReference:'legacy_nominal_reference.json',legacySourceSHA256:legacy.sourceSHA256,
 changes:['M3 forceMask controls force equations and active inverse parameters','Disabled coefficients become zero in theta and effective output config','M1 fit rejection reports actual event/domain reason'],
 scope:'Code and declared-equation verification, not evidence that the enabled force combination is physically true.'};
if(process.argv.includes('--write-report'))fs.writeFileSync(path.join(__dirname,'force-mask.validation.json'),JSON.stringify(report,null,2));
console.log(JSON.stringify({passed:true,checks:checks.length,version:P.VERSION,independentMagnusSyntheticMaximumParameterError:syntheticError},null,2));
