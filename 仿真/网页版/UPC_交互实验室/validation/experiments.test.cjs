'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const ROOT=fs.existsSync(path.resolve(__dirname,'../dist/physics.js'))?path.resolve(__dirname,'../dist'):path.resolve(__dirname,'..');
const {createPhysics}=require(path.join(ROOT,'physics.js')),{createDensity}=require(path.join(ROOT,'density.js')),{createStudy}=require(path.join(ROOT,'study.js')),{createMapping}=require(path.join(ROOT,'mapping.js')),{createExperiments}=require(path.join(ROOT,'experiments.js'));
const P=createPhysics(),D=createDensity(),S=createStudy(P,D),M=createMapping(P,D,S),E=createExperiments(P,D,S,M),checks=[];
const clone=x=>JSON.parse(JSON.stringify(x)),close=(a,b,t=1e-10)=>assert(Math.abs(a-b)<=t,a+' != '+b),check=(name,fn)=>{fn();checks.push({name,passed:true});};
const cfg=E.defaultConfig(),study={...S.defaultStudyConfig(),sampleCount:8,assessmentCount:8,intervalReplicates:50,schemes:S.generateSchemes({allowDrag:false,allowMagnus:false,includeConstantAccel:false}),selection:{...S.defaultStudyConfig().selection,minSamples:2}};
check('Forward reference generation, fine step and seed reproducibility',()=>{
 for(const model of ['M0','M2','M3']){const a=E.generateReference({...cfg.reference,model},.01,123),b=E.generateReference({...cfg.reference,model},.01,123);assert.equal(a.status,'ok');assert.deepEqual(a,b);assert.equal(a.referenceConfig.h,.0005);assert.equal(a.referenceConfig.fitH,.0005);assert.equal(a.provenance.independentReferenceAlgorithm,false);a.observations.forEach((o,i)=>o.position.forEach((x,j)=>close(x,a.provenance.noiseFreeObservations[i].position[j]+a.provenance.noise_m[i][j])));}
 const a=E.generateReference({...cfg.reference,model:'M0'},0,123);a.observations.forEach(o=>o.position.forEach((x,j)=>close(x,cfg.reference.initial[j]+cfg.reference.velocity[j]*o.t+(j===2?-.5*cfg.reference.g*o.t**2:0),1e-11)));
});
check('Only observations cross into the inverse configuration',()=>{
 let captured;const spy={...M,predictFromFrames:(a,b,c)=>{captured={a:clone(a),b:clone(b),c:clone(c)};return M.predictFromFrames(a,b,c);}},A=createExperiments(P,D,S,spy);
 const truth={...cfg.reference,m:.46,R:.1,rho:1.1,wind:[5,-3,.2],velocity:[32,4,7],Cd:.4,b:[0,.2,-.3]},r=A.runExperiment(truth,study,{observationSigma:.01,seed:246});
 assert.equal(r.status,'ok');assert.deepEqual(captured.c,study);assert.deepEqual(captured.a,r.reference.observations[0].position);assert.deepEqual(captured.b,r.reference.observations[1].position);assert.equal(r.study.config.base.m,study.base.m);assert.notEqual(r.study.config.base.m,truth.m);assert.equal(r.provenance.physicalValidation,false);
});
check('Manual edits change inverse results but not the reference trajectory',()=>{
 const a=E.runExperiment(cfg.reference,study,{observationSigma:.01,seed:987}),frames=clone(a.reference.observations);frames[0].position[1]+=.03;
 const b=E.runExperiment(cfg.reference,study,{observationSigma:.01,seed:987,observationsOverride:frames});assert.equal(b.status,'ok');assert.deepEqual(a.reference.trajectory,b.reference.trajectory);assert.deepEqual(a.reference.prediction,b.reference.prediction);assert.equal(b.provenance.manualObservationOverride,true);assert.deepEqual(b.provenance.generatedObservations,a.reference.observations);assert.deepEqual(b.study.config.base.observations,frames);assert.notDeepEqual(a.study.schemes[0].nominal.theta,b.study.schemes[0].nominal.theta);
});
check('Empirical ellipse mean, covariance, area and membership have an analytic reference',()=>{
 const points=Array.from({length:100},(_,i)=>({y:2*Math.cos(2*Math.PI*i/100),z:1+Math.sin(2*Math.PI*i/100)})),r=E.empiricalRegion(points,{y:0,z:1},.95,6),out=E.empiricalRegion(points,{y:3,z:1},.95,6);
 assert.equal(r.status,'ok');close(r.mean[0],0);close(r.mean[1],1);close(r.covariance[0][0],200/99);close(r.covariance[1][1],50/99);close(r.threshold,1.98,1e-12);close(r.area_m2,2*Math.PI,1e-10);assert(r.truthInside);assert.equal(out.truthInside,false);assert.equal(r.calibrationStatus,'not_calibrated');assert.equal(E.empiricalRegion(Array(8).fill({y:1,z:1}),null).status,'degenerate');assert.equal(E.empiricalRegion([{y:1,z:1}],null).status,'insufficient_samples');
});
check('Every candidate uses the same reference and labels cloud provenance',()=>{
 const c={...study,schemes:S.generateSchemes()},r=E.runExperiment(cfg.reference,c,{observationSigma:.01,seed:123});assert.equal(r.evaluations.length,c.schemes.length);
 r.evaluations.forEach(e=>{assert.equal(e.referenceStatus,r.reference.prediction.status);assert(['assessment','pilot'].includes(e.regionSource));assert.equal(e.cloudTotalTrials,8);if(e.estimatedArrival&&e.truthArrival){const p=e.nominalPrediction,t=r.reference.prediction;close(e.positionError_m,Math.hypot(p.yGoal-t.yGoal,p.zGoal-t.zGoal));close(e.timeError_s,p.tGoal-t.tGoal);}});assert.equal(r.calibration.status,'not_calibrated');
});
check('Batch failures retain the full denominator and omit bulky trajectories',()=>{
 const c={...cfg,study,repetitions:1,seed:456,scan:{axes:[{path:'velocity.0',values:[30,120]}]}},progress=[],r=E.runBatch(c,x=>progress.push(x));assert.equal(r.samples.length,2);assert.equal(r.samples[1].status,'reference_failed');assert.equal(r.matrix.caseIds.length,2);
 r.summary.forEach(s=>{assert.equal(s.totalCases,2);assert.equal(s.evaluatedCases,1);assert.equal(s.unavailableCases,1);assert.equal(s.predictionStatusCounts.not_run,1);});assert(r.samples.every(s=>!Object.hasOwn(s.reference,'trajectory')));assert.equal(progress.at(-1).completed,2);assert.deepEqual(r,E.runBatch(c));
});
check('Two axes, independent test seeds and unsupported calibration are explicit',()=>{
 const fast={...study,sampleCount:2,assessmentCount:2},c={...cfg,study:fast,repetitions:1,scan:{axes:[{path:'Cd',values:[.2,.4]},{path:'forceMask.magnus',values:[false,true]}]}};assert.equal(E.runBatch(c).samples.length,4);
 const one={...cfg,study:fast,repetitions:1,scan:{axes:[]}},a=E.runBatch(one),b=E.runBatch({...one,phase:'test'});assert.notEqual(a.samples[0].seed,b.samples[0].seed);assert.notDeepEqual(a.samples[0].reference.observations,b.samples[0].reference.observations);assert.equal(b.design.heldOutCalibration,false);assert.equal(b.calibration.status,'not_calibrated');assert.throws(()=>E.runBatch({...one,calibration:{enabled:true}}));assert.throws(()=>E.runBatch({...one,phase:'test',testSeed:one.seed}));
});
check('Zero noise and unavailable first-flight observations do not fabricate uncertainty',()=>{
 const r=E.runExperiment(cfg.reference,{...study,sigma:0},{observationSigma:0,seed:1});assert(r.evaluations.every(e=>e.region.status!=='ok'));const failed=E.generateReference({...cfg.reference,velocity:[30,5,-1]},.01,1);assert.equal(failed.status,'reference_observations_unavailable');assert.equal(failed.observations,null);
});
check('All factories serialize in a Worker without external dependencies',()=>{
 const code='(()=>{const P=('+createPhysics.toString()+')(),D=('+createDensity.toString()+')(),S=('+createStudy.toString()+')(P,D),M=('+createMapping.toString()+')(P,D,S);return ('+createExperiments.toString()+')(P,D,S,M);})()',W=vm.runInNewContext(code);assert.deepEqual(clone(W.generateReference(clone(cfg.reference),.01,123)),E.generateReference(cfg.reference,.01,123));
});
const report={passed:true,checks:checks.length,version:E.VERSION,node:process.version,tests:checks,scope:'Small software fixtures, not paper experiments; no calibrated coverage or real-ball validation.'};
if(process.argv.includes('--write-report'))fs.writeFileSync(path.join(__dirname,'experiments.validation.json'),JSON.stringify(report,null,2));console.log(JSON.stringify(report,null,2));
