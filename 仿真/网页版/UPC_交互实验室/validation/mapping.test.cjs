'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),crypto=require('node:crypto');
const {createPhysics}=require('../physics.js'),{createDensity}=require('../density.js'),{createStudy}=require('../study.js'),{createMapping}=require('../mapping.js');
const P=createPhysics(),D=createDensity(),S=createStudy(P,D),M=createMapping(P,D,S),tests=[];
const clone=x=>JSON.parse(JSON.stringify(x));
const equal=(a,b)=>assert.deepEqual(clone(a),clone(b));
const check=(name,f)=>{f();tests.push({name,passed:true});};
const config={...S.defaultStudyConfig(),sampleCount:20,assessmentCount:20,intervalReplicates:50,seed:13579};
const p05=[1.5,.45,.44],p10=[2.96,.86,.74],changed05=[1.5,.48,.44];
const first=M.predictFromFrames(p05,p10,config),changed=M.predictFromFrames(changed05,p10,config);

check('New two-frame coordinates rerun every inverse fit and produce new predictions and densities',()=>{
 equal(first.config.base.observations.map(o=>o.position),[p05,p10]);
 equal(changed.config.base.observations.map(o=>o.position),[changed05,p10]);
 assert.deepEqual(first.schemes.map(s=>s.id),['M0','M1','M2','M3','G_M']);
 changed.schemes.forEach((s,i)=>{
  const independent=P.solve({...s.physicsConfig,observations:[{t:.05,position:changed05},{t:.1,position:p10}]});
  equal(s.nominal.theta,independent.theta);equal(s.nominal.prediction,independent.prediction);
  assert.notDeepEqual(s.nominal.theta,first.schemes[i].nominal.theta);
  assert.notDeepEqual(s.nominal.prediction,first.schemes[i].nominal.prediction);
  if(s.density.status==='ok'&&first.schemes[i].density.status==='ok')assert.notDeepEqual(s.density.centers,first.schemes[i].density.centers);
 });
});
check('Mapped selection is linked to the fresh nominal solve, with ambiguity and independent assessment retained',()=>{
 for(const r of [first,changed]){
  assert.equal(r.mappedPrediction.unique,r.decision.uniqueSupported);
  assert.equal(r.mappedPrediction.provisional,!r.decision.uniqueSupported);
  equal(r.mappedPrediction.candidates.map(s=>s.id),r.schemes.filter(s=>r.decision.candidateIds.includes(s.id)).map(s=>s.id));
  if(r.decision.preferredId){
   const selected=r.schemes.find(s=>s.id===r.decision.preferredId);
   assert.equal(r.mappedPrediction.selected.id,selected.id);equal(r.mappedPrediction.selected.prediction,selected.nominal.prediction);
   const assessed=r.assessment.models.find(s=>s.id===selected.id);
   assert.equal(r.uncertainty.selectedIndependentAssessment.schemeId,selected.id);
   equal(r.uncertainty.selectedIndependentAssessment.summary,assessed.summary);
   assert.equal(assessed.summary.total,config.assessmentCount);
  }
 }
});
check('Common perturbations are regenerated deterministically around the changed observations',()=>{
 first.trials.forEach((t,i)=>{
  equal(t.noise_m,changed.trials[i].noise_m);
  for(let j=0;j<2;j++)for(let k=0;k<3;k++){
   const expected=(j===0?changed05:p10)[k]+t.noise_m[3*j+k];
   assert(Math.abs(changed.trials[i].observations[j].position[k]-expected)<1e-14);
  }
 });
 assert.notDeepEqual(first.assessment.trials[0].noise_m,first.trials[0].noise_m);
});
check('Invalid or incompatible new frames never fall back to a previous saved density',()=>{
 assert.throws(()=>M.predictFromFrames([NaN,0,1],p10,config));
 const r=M.predictFromFrames([20,.45,.44],[40,.86,.74],config);
 assert.equal(r.mappedPrediction.selected,null);assert.equal(r.decision.preferredId,null);
 assert.equal(r.uncertainty.selectedIndependentAssessment,null);
 r.schemes.forEach(s=>assert.equal(s.summary.eligible,false));
});
let portableBytes=0;
check('CommonJS exportSource produces a standalone module in an empty VM with no files or globals',()=>{
 const source=M.exportSource(config);portableBytes=Buffer.byteLength(source);
 const context={module:{exports:{}}};vm.createContext(context);assert.equal(context.require,undefined);
 vm.runInContext(source,context);const standalone=context.module.exports;
 assert.equal(typeof standalone.F,'function');assert.equal(standalone.F,standalone.predictFromFrames);
 const original=standalone.F(p05,p10),modified=standalone.F(changed05,p10);
 equal(original.mappedPrediction,first.mappedPrediction);equal(original.decision,first.decision);
 equal(original.trials,first.trials);equal(original.assessment.models.map(s=>s.summary),first.assessment.models.map(s=>s.summary));
 equal(modified.mappedPrediction,changed.mappedPrediction);equal(modified.schemes.map(s=>s.density),changed.schemes.map(s=>s.density));
 assert.notDeepEqual(clone(original.schemes[0].nominal.theta),clone(modified.schemes[0].nominal.theta));
});
check('Delivered single-file mapping artifact also executes fresh observations without external imports',()=>{
 const filename=path.join(__dirname,'../综合预测映射.js');assert(fs.existsSync(filename));
 const context={module:{exports:{}}};vm.createContext(context);vm.runInContext(fs.readFileSync(filename,'utf8'),context);
 const exported=context.module.exports,r=exported.F(changed05,p10,config);
 equal(r.mappedPrediction,changed.mappedPrediction);equal(r.decision,changed.decision);
});
let outer,outerCalls=[];
check('Outer selection sensitivity reruns the whole selector and avoids an extra assessment/noise layer',()=>{
 const watched={...S,runStudy:c=>{outerCalls.push(clone(c));return S.runStudy(c);}};
 const W=createMapping(P,D,watched);
 outer=W.selectionSensitivity(p05,p10,config,{outerSampleCount:3,seed:24680});
 assert.equal(outerCalls.length,3);assert.equal(outer.rows.length,3);assert.equal(outer.outerSampleCount,3);
 outerCalls.forEach((c,i)=>{assert.equal(c.assessmentEnabled,false);assert.equal(c.seed,config.seed);equal(c.base.observations,outer.rows[i].observations);});
 outer.rows.forEach(row=>{
  if(!row.selectedId)return;
  const scheme=config.schemes.find(s=>s.id===row.selectedId);
  const direct=P.solve({...S.schemeConfig(config,scheme),observations:row.observations});
  equal(row.prediction,direct.prediction);
 });
 const chosen=Object.values(outer.selectionCounts).reduce((a,b)=>a+b,0);
 assert.equal(chosen+outer.noEligibleCount,3);assert.equal(outer.arrivals,outer.rows.filter(row=>row.status==='crossed_goal_plane').length);
 assert.equal(outer.responseDensity.planeArrivalMass.nTotal,3);assert(outer.ambiguousCount<=chosen);
});
check('Outer seeds and all failure denominators reproduce without implying model probabilities',()=>{
 const again=M.selectionSensitivity(p05,p10,config,{outerSampleCount:3,seed:24680});equal(again,outer);
 const failed=M.selectionSensitivity([20,.45,.44],[40,.86,.74],config,{outerSampleCount:2,seed:123});
 assert.equal(failed.rows.length,2);assert.equal(failed.noEligibleCount,2);assert.equal(failed.arrivals,0);
 assert.equal(failed.responseDensity.planeArrivalMass.nTotal,2);assert.equal(failed.responseDensity.status,'degenerate');
 assert(Object.values(failed.selectionCounts).every(n=>n===0));
 assert.throws(()=>M.selectionSensitivity(p05,p10,config,{outerSampleCount:1}));
});
const sources=Object.fromEntries(['physics.js','density.js','study.js','mapping.js','综合预测映射.js'].map(name=>['../'+name,crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'..',name))).digest('hex')]));
const report={passed:true,checks:tests.length,node:process.version,mappingVersion:M.VERSION,tests,sources,
 sampleCount:config.sampleCount,assessmentCount:config.assessmentCount,seed:config.seed,portableBytes,
 outerTest:{sampleCount:3,seed:24680,selectionCounts:outer.selectionCounts,ambiguousCount:outer.ambiguousCount,noEligibleCount:outer.noEligibleCount},
 scope:'Software mapping and denominator checks only. The tiny outer test budget is not a research estimate or evidence of real physical model probabilities.'};
if(process.argv.includes('--write-report'))fs.writeFileSync(path.join(__dirname,'mapping.validation.json'),JSON.stringify(report,null,2));
console.log(JSON.stringify({passed:true,checks:tests.length,portableBytes,outerCounts:outer.selectionCounts},null,2));
