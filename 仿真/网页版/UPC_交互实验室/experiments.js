/* Simulated-reference laboratory. Serializable factory; no DOM or external imports. */
function createExperiments(P,D,S,M){
'use strict';
const VERSION='upc-experiments-0.1',clone=x=>JSON.parse(JSON.stringify(x)),finite3=x=>Array.isArray(x)&&x.length===3&&x.every(Number.isFinite);
const arrived=p=>p?.status==='crossed_goal_plane'&&[p.yGoal,p.zGoal,p.tGoal].every(Number.isFinite);
function random(seed){let a=seed>>>0;return()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function gaussian(r){return Math.sqrt(-2*Math.log(Math.max(1e-15,r())))*Math.cos(2*Math.PI*r());}
function quantile(v,p){const a=v.slice().sort((x,y)=>x-y),u=p*(a.length-1),i=Math.floor(u);return a[i]+(a[Math.min(i+1,a.length-1)]-a[i])*(u-i);}
function defaultConfig(){return {version:VERSION,reference:{...P.defaultConfig(),mode:'forward',model:'M3',velocity:[30,5,6],Cd:.3,b:[0,.1,-.25],forceMask:{drag:true,magnus:true},h:.0005,fitH:.0005},study:S.defaultStudyConfig(),observationSigma:.01,seed:20261107,testSeed:20261108,phase:'exploration',repetitions:3,scan:{axes:[]},regionLevel:.95,regionMinimumSamples:6,calibration:{enabled:false,status:'not_calibrated'}};}
function calibration(level){return {status:'not_calibrated',nominalLevel:level,interpretation:'Empirical sample-cloud ellipse; no trained/calibrated coverage guarantee. An independently seeded test batch is not automatically a calibrated test.'};}
function generateReference(referenceConfig,observationSigma=.01,seed=20261107){
 if(!referenceConfig||typeof referenceConfig!=='object'||Array.isArray(referenceConfig))throw Error('referenceConfig must be an object.');
 if(!Number.isFinite(observationSigma)||observationSigma<0||!Number.isSafeInteger(seed))throw Error('Noise sigma must be nonnegative; seed must be a safe integer.');
 const cfg={...P.defaultConfig(),...clone(referenceConfig),mode:'forward',h:.0005,fitH:.0005};
 cfg.observations=[.05,.1].map(t=>({t,position:finite3(cfg.initial)?cfg.initial.slice():[0,0,.11]}));
 const provenance={kind:'simulated_reference',physicsVersion:P.VERSION||null,referenceStep_s:.0005,independentReferenceAlgorithm:false,physicalValidation:false,observationSigma_m:observationSigma,noiseSeed:seed,manualObservationOverride:false,noiseModel:'Six independent additive Gaussian coordinate errors; no clipping.',numericalMethod:'Same forward RK4 kernel at 0.0005 s with event bisection; finer simulated reference, not independent algorithm or real-ball truth.',hiddenParameterPolicy:'Only observed positions enter the inverse mapping; reference parameters never overwrite declared study inputs.'};
 let run;try{run=P.solve(cfg);}catch(error){return {status:'reference_failed',trajectory:[],prediction:{status:'integration_failed'},observations:null,referenceConfig:cfg,provenance,issues:[String(error)]};}
 const result={status:'ok',trajectory:run.trajectory||[],prediction:run.prediction,observations:null,referenceConfig:cfg,provenance,issues:run.issues||[]};
 if(!run.prediction||['invalid_input','fit_failed','integration_failed','outside_model_domain'].includes(run.prediction.status))return {...result,status:'reference_failed'};
 if((run.prediction.status==='ground_before_goal'&&run.prediction.eventTime<.1-1e-10)||(arrived(run.prediction)&&run.prediction.tGoal<=.1))return {...result,status:'reference_observations_unavailable',issues:[...result.issues,'First-flight contact or target event precedes the last observation.']};
 let clean;try{clean=P.predictPositions(cfg,run.theta,[.05,.1]);}catch(error){return {...result,status:'reference_failed',issues:[...result.issues,String(error)]};}
 if(clean.some(p=>!finite3(p)))return {...result,status:'reference_failed',issues:[...result.issues,'Nonfinite reference observation.']};
 const rng=random(seed),noise=clean.map(()=>[observationSigma*gaussian(rng),observationSigma*gaussian(rng),observationSigma*gaussian(rng)]);
 const truthFrames=clean.map((position,i)=>({t:[.05,.1][i],position:position.slice()}));
 return {...result,observations:clean.map((p,i)=>({t:[.05,.1][i],position:p.map((v,j)=>v+noise[i][j])})),referenceConfig:{...cfg,observations:truthFrames},provenance:{...provenance,noise_m:noise,noiseFreeObservations:truthFrames}};
}
function empiricalRegion(points,truth,level=.95,minimumSamples=6){
 const base={status:'insufficient_samples',nominalLevel:level,sampleCount:points.length,mean:null,covariance:null,threshold:null,area_m2:null,truthInside:null,truthMahalanobisSquared:null,calibrationStatus:'not_calibrated',geometry:'Full ellipse, not clipped at the ground or frame.',interpretation:'Empirical quantile of in-cloud squared Mahalanobis distances; NOT a calibrated future-coverage region.'};
 if(points.length<minimumSamples)return base;
 if(points.some(p=>!Number.isFinite(p.y)||!Number.isFinite(p.z)))return {...base,status:'invalid_samples'};
 const n=points.length,mean=[points.reduce((a,p)=>a+p.y,0)/n,points.reduce((a,p)=>a+p.z,0)/n];let yy=0,yz=0,zz=0;
 for(const p of points){const y=p.y-mean[0],z=p.z-mean[1];yy+=y*y;yz+=y*z;zz+=z*z;}yy/=n-1;yz/=n-1;zz/=n-1;
 const det=yy*zz-yz*yz,trace=yy+zz,partial={...base,mean,covariance:[[yy,yz],[yz,zz]]};
 if(!Number.isFinite(det)||trace<=0||det<=1e-12*trace*trace)return {...partial,status:'degenerate'};
 const score=p=>{const y=p.y-mean[0],z=p.z-mean[1];return Math.max(0,(zz*y*y-2*yz*y*z+yy*z*z)/det);};
 const scores=points.map(score),threshold=quantile(scores,level),area=Math.PI*threshold*Math.sqrt(det),truthScore=truth&&Number.isFinite(truth.y)&&Number.isFinite(truth.z)?score(truth):null;
 if(!Number.isFinite(area)||area<=0)return {...partial,status:'degenerate'};
 return {...partial,status:'ok',threshold,area_m2:area,empiricalIncludedCount:scores.filter(x=>x<=threshold).length,truthInside:truthScore===null?null:truthScore<=threshold,truthMahalanobisSquared:truthScore};
}
function optionsRegion(o){const level=o.regionLevel??.95,min=o.regionMinimumSamples??6;if(!Number.isFinite(level)||level<=0||level>=1||!Number.isSafeInteger(min)||min<3)throw Error('Region level must be in (0,1), minimum samples an integer >=3.');return {level,min};}
function unavailable(s,reference,reason,level){return {schemeId:s.id,label:s.label||s.id,status:'unavailable',reason,referenceStatus:reference.prediction?.status||'unavailable',predictionStatus:'not_run',eventMatch:null,truthArrival:arrived(reference.prediction),estimatedArrival:false,positionError_m:null,timeError_s:null,absoluteTimeError_s:null,regionSource:'none',nominalPrediction:null,region:empiricalRegion([],null,level)};}
function evaluate(s,study,reference,level,min){
 const assessed=study.assessment?.models?.find(a=>a.id===s.id),cloud=assessed||s,samples=cloud.samples||[],points=samples.filter(p=>p.status==='crossed_goal_plane').map(p=>({y:p.yGoal,z:p.zGoal}));
 const truth=reference.prediction,p=s.nominal.prediction,both=arrived(truth)&&arrived(p);
 return {schemeId:s.id,label:s.label||s.id,status:'evaluated',referenceStatus:truth.status,predictionStatus:p.status,eventMatch:truth.status===p.status,truthArrival:arrived(truth),estimatedArrival:arrived(p),positionError_m:both?Math.hypot(p.yGoal-truth.yGoal,p.zGoal-truth.zGoal):null,timeError_s:both?p.tGoal-truth.tGoal:null,absoluteTimeError_s:both?Math.abs(p.tGoal-truth.tGoal):null,nominalPrediction:clone(p),regionSource:assessed?'assessment':'pilot',regionPurpose:assessed?'Separate perturbation stream after selection, not trajectory-level calibration.':'Reuses selection samples; descriptive and potentially optimistic.',cloudTotalTrials:samples.length,cloudArrivalCount:points.length,cloudStatusCounts:clone(cloud.summary?.statusCounts||{}),region:empiricalRegion(points,arrived(truth)?{y:truth.yGoal,z:truth.zGoal}:null,level,min)};
}
function runExperiment(referenceConfig,declaredStudyConfig,options={},onProgress){
 const {level,min}=optionsRegion(options),studyConfig=clone(declaredStudyConfig||S.defaultStudyConfig()),sigma=options.observationSigma??studyConfig.sigma??.01,seed=options.seed??20261107;
 const reference=generateReference(referenceConfig,sigma,seed),schemes=studyConfig.schemes||S.defaultStudyConfig().schemes;
 if(onProgress)onProgress({stage:'reference',completed:1,total:2,message:'模拟参考已生成，反演仅接收两帧位置。'});
 if(reference.status!=='ok')return {version:VERSION,status:'reference_failed',reference,study:null,declaredStudyConfig:studyConfig,evaluations:schemes.map(s=>unavailable(s,reference,reference.status,level)),provenance:reference.provenance,calibration:calibration(level)};
 if(options.observationsOverride!==undefined&&options.observationsOverride!==null){
  const o=options.observationsOverride;if(!Array.isArray(o)||o.length!==2||o.some((p,i)=>!p||p.t!==[.05,.1][i]||!finite3(p.position)))throw Error('observationsOverride requires two finite {t,position} frames at 0.05 and 0.10 seconds.');
  reference.provenance.generatedObservations=clone(reference.observations);reference.observations=clone(o);reference.provenance.manualObservationOverride=true;reference.provenance.observationInterpretation='Manually edited observation scenario; the simulated reference trajectory and parameters are unchanged. Gaussian coverage interpretation does not apply to manual edits.';
 }
 let study;try{
  study=M.predictFromFrames(reference.observations[0].position,reference.observations[1].position,studyConfig);
 }catch(error){return {version:VERSION,status:'study_failed',reference,study:null,declaredStudyConfig:studyConfig,evaluations:schemes.map(s=>unavailable(s,reference,String(error),level)),issues:[String(error)],provenance:reference.provenance,calibration:calibration(level)};}
 const evaluations=study.schemes.map(s=>evaluate(s,study,reference,level,min));
 if(onProgress)onProgress({stage:'evaluation',completed:2,total:2,message:'已计算每个候选的参考误差与未校准经验区域。'});
 return {version:VERSION,status:'ok',reference,study,declaredStudyConfig:studyConfig,evaluations,provenance:reference.provenance,calibration:calibration(level)};
}
const scanPaths=['model','velocity.0','velocity.1','velocity.2','Cd','b.0','b.1','b.2','wind.0','wind.1','wind.2','m','R','rho','g','initial.0','initial.1','initial.2','constantAcceleration.0','constantAcceleration.1','constantAcceleration.2','forceMask.drag','forceMask.magnus'];
function setPath(o,path,value){const keys=path.split('.');for(let i=0;i<keys.length-1;i++)o=o[keys[i]];o[keys[keys.length-1]]=value;}
function batchConfig(input){
 if(input!==undefined&&(!input||typeof input!=='object'||Array.isArray(input)))throw Error('Batch config must be an object.');
 const d=defaultConfig(),c={...d,...clone(input||{})};c.reference={...d.reference,...c.reference,forceMask:{...d.reference.forceMask,...c.reference?.forceMask}};c.study={...d.study,...c.study};
 const axes=c.scan?.axes||[];if(!Array.isArray(axes)||axes.length>2||new Set(axes.map(a=>a.path)).size!==axes.length)throw Error('Use up to two distinct scan axes.');
 for(const a of axes){if(!scanPaths.includes(a.path)||!Array.isArray(a.values)||!a.values.length||a.values.length>100)throw Error('Invalid scan path or values.');
  if(a.values.some(v=>a.path==='model'?!['M0','M1','M2','M3'].includes(v):a.path.startsWith('forceMask.')?typeof v!=='boolean':!Number.isFinite(v)))throw Error('Invalid scan value type.');}
 if(!Number.isSafeInteger(c.repetitions)||c.repetitions<1||!Number.isSafeInteger(c.seed)||!Number.isSafeInteger(c.testSeed))throw Error('Positive integer repetitions and safe integer seeds required.');
 if(!['exploration','test'].includes(c.phase))throw Error('phase must be exploration or test; calibration is not implemented.');
 if(c.phase==='test'&&c.testSeed===c.seed)throw Error('testSeed must differ from exploration seed.');
 if(!Number.isFinite(c.observationSigma)||c.observationSigma<0)throw Error('Observation sigma must be nonnegative.');
 const caseCount=axes.reduce((n,a)=>n*a.values.length,c.repetitions);if(caseCount>500)throw Error('Interactive compute budget is at most 500 cases.');
 if(c.calibration?.enabled)throw Error('Coverage calibration is not implemented; no calibrated report can be requested.');
 optionsRegion(c);return {...c,scan:{axes},caseCount};
}
function summarize(samples,ids){return ids.map(schemeId=>{
 const rows=samples.map(s=>s.evaluations.find(e=>e.schemeId===schemeId)),evaluated=rows.filter(r=>r?.status==='evaluated'),truth=rows.filter(r=>r?.truthArrival),paired=rows.filter(r=>Number.isFinite(r?.positionError_m)),regions=truth.filter(r=>r.region.status==='ok'),covered=regions.filter(r=>r.region.truthInside===true),counts={};
 for(const r of rows){const key=r?.predictionStatus||'missing';counts[key]=(counts[key]||0)+1;}
 return {schemeId,totalCases:samples.length,evaluatedCases:evaluated.length,unavailableCases:samples.length-evaluated.length,truthArrivalCases:truth.length,bothArrivalCases:paired.length,regionAvailableTruthArrivalCases:regions.length,coveredTruthArrivalCases:covered.length,coverageAmongAvailableRegions:regions.length?covered.length/regions.length:null,coverageAmongAllTruthArrivals:truth.length?covered.length/truth.length:null,coveredFractionAllAttemptedCases:covered.length/samples.length,positionRMSE_m:paired.length?Math.sqrt(paired.reduce((s,r)=>s+r.positionError_m**2,0)/paired.length):null,meanAbsoluteTimeError_s:paired.length?paired.reduce((s,r)=>s+r.absoluteTimeError_s,0)/paired.length:null,meanRegionArea_m2:regions.length?regions.reduce((s,r)=>s+r.region.area_m2,0)/regions.length:null,eventMatchFractionAmongEvaluated:evaluated.length?evaluated.filter(r=>r.eventMatch).length/evaluated.length:null,predictionStatusCounts:counts,regionSources:{assessment:rows.filter(r=>r?.regionSource==='assessment').length,pilot:rows.filter(r=>r?.regionSource==='pilot').length},calibrationStatus:'not_calibrated',denominatorNote:'All-truth-arrival coverage counts unavailable regions as misses; all-attempted fraction also includes nonarrival/reference failures and is not spatial coverage.'};
});}
function runBatch(input,onProgress){
 const config=batchConfig(input),combinations=[{}];for(const axis of config.scan.axes){const old=combinations.splice(0);for(const p of old)for(const value of axis.values)combinations.push({...p,[axis.path]:value});}
 const samples=[],schemeIds=config.study.schemes.map(s=>s.id),phaseSeed=config.phase==='test'?config.testSeed:config.seed;let i=0;
 for(const parameters of combinations)for(let repetition=0;repetition<config.repetitions;repetition++){
  const ref=clone(config.reference);for(const [key,value]of Object.entries(parameters))setPath(ref,key,value);
  const seed=(phaseSeed+Math.imul(i+1,0x9e3779b1))>>>0;
  const r=runExperiment(ref,config.study,{observationSigma:config.observationSigma,seed,regionLevel:config.regionLevel,regionMinimumSamples:config.regionMinimumSamples});
  samples.push({index:i,caseId:'case_'+String(i+1).padStart(4,'0'),phase:config.phase,repetition,parameters,seed,status:r.status,reference:{status:r.reference.status,prediction:r.reference.prediction,observations:r.reference.observations,referenceConfig:r.reference.referenceConfig,provenance:r.reference.provenance,issues:r.reference.issues},decision:r.study?.decision||null,evaluations:r.evaluations,issues:r.issues||[]});i++;
  if(onProgress)onProgress({stage:'batch',completed:i,total:config.caseCount,caseIndex:i,caseCount:config.caseCount,message:'模拟参考验证 '+i+'/'+config.caseCount});
 }
 const cell=(s,id)=>s.evaluations.find(e=>e.schemeId===id),matrix={caseIds:samples.map(s=>s.caseId),schemeIds,layout:'rows=cases; columns=candidates'};
 for(const field of ['positionError_m','eventMatch','predictionStatus'])matrix[field]=samples.map(s=>schemeIds.map(id=>cell(s,id)?.[field]??null));
 for(const [field,key]of [['area_m2','area_m2'],['truthInside','truthInside']])matrix[field]=samples.map(s=>schemeIds.map(id=>cell(s,id)?.region[key]??null));
 return {version:VERSION,status:'completed',config,samples,matrix,summary:summarize(samples,schemeIds),calibration:calibration(config.regionLevel),design:{pairedCases:true,phase:config.phase,phaseSeed,testSeed:config.testSeed,referenceStep_s:.0005,independentReferenceAlgorithm:false,heldOutCalibration:false,seedPolicy:'Distinct deterministic observation seed per case; declared inner study seed fixed. Test phase uses separate testSeed, but candidate/tuning freeze is the operator responsibility.',scope:'Simulated-reference comparisons only; no empirical population or real-ball validation.'}};
}
return {VERSION,defaultConfig,generateReference,runExperiment,runBatch,empiricalRegion};
}
if(typeof module!=='undefined'&&module.exports)module.exports={createExperiments};
if(typeof globalThis!=='undefined')globalThis.createExperiments=createExperiments;

