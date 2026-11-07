/* A portable input-to-prediction mapping, including candidate enumeration and uncertainty.
 * No fixed fitted density is substituted for a fresh inverse calculation.
 */
function createMapping(P,D,S){
 'use strict';
 const VERSION='upc-prediction-mapping-0.2',clone=x=>JSON.parse(JSON.stringify(x));
 function configFor(p05,p10,input){
  for(const p of [p05,p10])if(!Array.isArray(p)||p.length!==3||!p.every(Number.isFinite))throw Error('两帧位置须分别为三个有限坐标，单位m');
  const base=S.defaultStudyConfig(),given=clone(input||{}),c={...base,...given};
  c.base={...base.base,...given.base,observations:[{t:.05,position:p05.slice()},{t:.1,position:p10.slice()}]};
  return c;
 }
 function compact(s){return {id:s.id,label:s.label,prediction:clone(s.nominal.prediction),theta:clone(s.nominal.theta),constantAcceleration:s.nominal.fit?.constantAcceleration??null};}
 function withMapping(study){
  const selected=study.schemes.find(s=>s.id===study.decision.preferredId),candidate=study.schemes.filter(s=>study.decision.candidateIds.includes(s.id));
  const selectedAssessment=study.assessment.models.find(s=>s.id===study.decision.preferredId);
  return {...study,mappingVersion:VERSION,mappedPrediction:{input:study.config.base.observations,selected:selected?compact(selected):null,candidates:candidate.map(compact),provisional:!study.decision.uniqueSupported,unique:study.decision.uniqueSupported,
   meaning:'选定组合在未扰动两帧观测上的首飞球心过门平面预测；候选并列时selected仅为暂定数值最小者。'},
   uncertainty:{selectedIndependentAssessment:selectedAssessment?{schemeId:selectedAssessment.id,summary:selectedAssessment.summary,density:selectedAssessment.density}:null,candidateAssessments:study.assessment.models.map(s=>({id:s.id,summary:s.summary})),
    betweenCandidateNominal:S.moments(candidate.filter(s=>s.nominal.prediction.status==='crossed_goal_plane').map(s=>({y:s.nominal.prediction.yGoal,z:s.nominal.prediction.zGoal}))),
    interpretation:'独立评估给定选定组合的输入扰动传播；候选间差异另列，不把组合枚举当实测概率。包含重新选型的敏感性需调用selectionSensitivity。'}};
 }
 function predictFromFrames(p05,p10,config){return withMapping(S.runStudy(configFor(p05,p10,config)));}
 function selectionSensitivity(p05,p10,config,options={}){
  const c=configFor(p05,p10,config),n=options.outerSampleCount??20;
  if(!Number.isInteger(n)||n<2||n>200)throw Error('外层重新选型次数须为2—200');
  const outerSeed=options.seed??((c.seed^0x7f4a7c15)>>>0);if(!Number.isInteger(outerSeed))throw Error('外层seed须为整数');
  const trials=S.makeTrials({...c,sampleCount:n,seed:outerSeed}),counts=Object.fromEntries(c.schemes.map(s=>[s.id,0])),rows=[],points=[];
  let ambiguous=0,noEligible=0;
  for(const trial of trials){
   const modified=configFor(trial.observations[0].position,trial.observations[1].position,c);
   // Fixed inner seed across outer inputs makes this a reproducible response of
   // the whole selector; changing only inner seeds would measure MC error instead.
   modified.assessmentEnabled=false;
   let r;try{r=S.runStudy(modified);}catch(error){rows.push({index:trial.index,status:'invalid_input',error:String(error),observations:trial.observations});noEligible++;continue;}
   const chosen=r.schemes.find(s=>s.id===r.decision.preferredId);
   if(!chosen){noEligible++;rows.push({index:trial.index,status:'no_eligible_scheme',decision:r.decision,observations:trial.observations});continue;}
   counts[chosen.id]++;if(!r.decision.uniqueSupported)ambiguous++;
   const pred=chosen.nominal.prediction;
   rows.push({index:trial.index,observations:trial.observations,selectedId:chosen.id,decision:r.decision,status:pred.status,prediction:pred});
   if(pred.status==='crossed_goal_plane')points.push({y:pred.yGoal,z:pred.zGoal});
  }
  return {version:VERSION,config:c,outerSampleCount:n,outerSeed,innerSeed:c.seed,rows,selectionCounts:counts,ambiguousCount:ambiguous,noEligibleCount:noEligible,arrivals:points.length,
   responseMoments:S.moments(points),responseDensity:D.buildDensity(points,{...c.bandwidth,supportZMin:c.base.R,planeArrivalMass:{nArrival:points.length,nTotal:n}}),
   interpretation:'外层扰动观测，每次重做完整候选比较，并仅保存该次选定模型在外层观测上的名义预测；无第二次叠加输出噪声。选中/并列频率是选择流程敏感性，不是真实模型概率，不是已校准真实覆盖。有限样本无法覆盖未纳入的物理机制。'};
 }
 function exportSource(config){
  const fixed=clone(config||S.defaultStudyConfig());
  // P, D and S factories are supplied when assembling this standalone module.
  const PF=typeof createPhysics==='function'?createPhysics:(typeof require==='function'?require('./physics.js').createPhysics:null),DF=typeof createDensity==='function'?createDensity:(typeof require==='function'?require('./density.js').createDensity:null),SF=typeof createStudy==='function'?createStudy:(typeof require==='function'?require('./study.js').createStudy:null);
  if(!PF||!DF||!SF)throw Error('导出需要可序列化的三个科学工厂');
  return `'use strict';\n// SI units. F([x05,y05,z05],[x10,y10,z10]) reruns inverse fits and model selection.\nconst createPhysics=${PF.toString()};\nconst createDensity=${DF.toString()};\nconst createStudy=${SF.toString()};\nconst createMapping=${createMapping.toString()};\nconst config=${JSON.stringify(fixed,null,2)};\nconst P=createPhysics(),D=createDensity(),S=createStudy(P,D),M=createMapping(P,D,S);\nfunction predictFromFrames(p05,p10,patch={}){return M.predictFromFrames(p05,p10,{...config,...patch,base:{...config.base,...patch.base}});}\nfunction selectionSensitivity(p05,p10,options={}){return M.selectionSensitivity(p05,p10,config,options);}\nconst F=predictFromFrames;\nconst exportedMapping={F,predictFromFrames,selectionSensitivity,config,version:M.VERSION};\nif(typeof module!=='undefined'&&module.exports)module.exports=exportedMapping;\nif(typeof globalThis!=='undefined')globalThis.UPCPredictionMapping=exportedMapping;\n`;
 }
 return {VERSION,predictFromFrames,selectionSensitivity,exportSource,withMapping};
}
if(typeof module!=='undefined'&&module.exports)module.exports={createMapping};
if(typeof globalThis!=='undefined')globalThis.createMapping=createMapping;
