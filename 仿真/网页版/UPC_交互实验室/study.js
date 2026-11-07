/* Shared-observation Monte Carlo comparison. No physics model is selected as physical truth. */
function createStudy(P,D){
 'use strict';
 const VERSION='upc-combination-study-0.2',clone=x=>JSON.parse(JSON.stringify(x));
 const STATES=['crossed_goal_plane','ground_before_goal','no_crossing_within_horizon','invalid_input','fit_failed','integration_failed','outside_model_domain'];
 const PRESETS={gravity:{model:'M0',label:'重力'},drag:{model:'M2',label:'重力＋阻力'},magnus:{model:'M3',label:'重力＋阻力＋Magnus'},wind:{model:'M3',label:'前三者＋给定风'},fields:{model:'M3',label:'受力组合'},constantAccel:{model:'M1',label:'常加速度代数对照'}};
 function generateSchemes(options={}){
  const {includeWind=false,includeConstantAccel=true,allowDrag=true,allowMagnus=true,wind=[0,2,0]}=options;
  const make=(id,preset,label,drag,magnus,usesWind=false)=>({id,preset,label,forceMask:{drag,magnus},usesWind});
  const out=[make('M0','gravity','M0 · 重力',false,false)];
  if(includeConstantAccel)out.push(make('M1','constantAccel','M1 · 常加速度代数对照',false,false));
  if(allowDrag)out.push(make('M2','drag','M2 · 重力＋阻力',true,false));
  if(allowDrag&&allowMagnus)out.push(make('M3','magnus','M3 · 重力＋阻力＋Magnus',true,true));
  if(allowMagnus)out.push(make('G_M','fields','重力＋Magnus（无阻力）',false,true));
  if(includeWind&&Array.isArray(wind)&&wind.some(v=>v!==0))for(const s of out.slice())if(s.forceMask.drag||s.forceMask.magnus)out.push({...clone(s),id:s.id+'_W',label:s.label+'＋指定风',usesWind:true});
  return out;
 }
 function defaultStudyConfig(){return {version:VERSION,base:P.defaultConfig(),sigma:.01,sampleCount:60,assessmentCount:60,assessmentEnabled:true,seed:20261107,wind:[0,2,0],schemes:generateSchemes(),bandwidth:{sigma_h:null,minBandwidth:.001},selection:{maxNormalizedRMSE:2,minArrivalFraction:.9,maxNumericalFailureFraction:.05,minSamples:20},intervalReplicates:200};}
 function normalize(input){const d=defaultStudyConfig(),c={...d,...clone(input||{})};c.base={...d.base,...(c.base||{})};c.bandwidth={...d.bandwidth,...(c.bandwidth||{})};c.selection={...d.selection,...(c.selection||{})};return c;}
 function validate(c){
  const issues=[];
  if(!Number.isFinite(c.sigma)||c.sigma<0||c.sigma>.5)issues.push('sigma须为0—0.5m的有限假设值');
  if(!Number.isInteger(c.sampleCount)||c.sampleCount<2||c.sampleCount>1000)issues.push('样本数须为2—1000整数');
  if(!Number.isInteger(c.assessmentCount)||c.assessmentCount<2||c.assessmentCount>1000)issues.push('独立评估样本数须为2—1000整数');
  if(typeof c.assessmentEnabled!=='boolean')issues.push('assessmentEnabled须为布尔值');
  if(!Number.isInteger(c.seed))issues.push('随机种子须为整数');
  if(!Array.isArray(c.wind)||c.wind.length!==3||!c.wind.every(Number.isFinite))issues.push('风须为3个有限数值');
  if(!Array.isArray(c.schemes)||c.schemes.length<1||c.schemes.length>32||c.schemes.some(s=>!s||!PRESETS[s.preset]||typeof s.id!=='string'||s.forceMask&&(!['drag','magnus'].every(k=>typeof s.forceMask[k]==='boolean')))||new Set(c.schemes.map(s=>s?.id)).size!==c.schemes.length)issues.push('须定义1—32个不同ID的可用候选组合，forceMask包含布尔drag/magnus');
  if(!Number.isFinite(c.selection.maxNormalizedRMSE)||c.selection.maxNormalizedRMSE<=0)issues.push('残差门限须为正数');
  for(const k of ['minArrivalFraction','maxNumericalFailureFraction'])if(!Number.isFinite(c.selection[k])||c.selection[k]<0||c.selection[k]>1)issues.push(k+'须在0—1之间');
  if(!Number.isInteger(c.selection.minSamples)||c.selection.minSamples<2)issues.push('最少样本门限须不小于2');
  if(!Number.isInteger(c.intervalReplicates)||c.intervalReplicates<50||c.intervalReplicates>2000)issues.push('区间重抽样次数须为50—2000');
  if(c.bandwidth.sigma_h!==null&&(!Number.isFinite(c.bandwidth.sigma_h)||c.bandwidth.sigma_h<=0))issues.push('平滑带宽须为空或正数');
  if(!Number.isFinite(c.bandwidth.minBandwidth)||c.bandwidth.minBandwidth<=0)issues.push('带宽下限须为正有限数');
  issues.push(...P.validateConfig({...c.base,model:'M3',mode:'inverse',sigma:c.sigma}));
  return issues;
 }
 function rng(seed){let a=seed>>>0;return()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return ((t^t>>>14)>>>0)/4294967296;};}
 function normal(random){return Math.sqrt(-2*Math.log(Math.max(1e-15,random())))*Math.cos(2*Math.PI*random());}
 function makeTrials(c){return Array.from({length:c.sampleCount},(_,index)=>{const seed=(c.seed+Math.imul(index+1,0x9e3779b1))>>>0,r=rng(seed),noise=Array.from({length:6},()=>c.sigma*normal(r));return {index,seed,noise_m:noise,observations:c.base.observations.map((o,j)=>({t:o.t,position:o.position.map((x,k)=>x+noise[3*j+k])}))};});}
 function schemeConfig(c,scheme){return {...clone(c.base),mode:'inverse',model:PRESETS[scheme.preset].model,sigma:c.sigma,forceMask:scheme.forceMask||{drag:true,magnus:true},wind:scheme.preset==='wind'||scheme.usesWind?c.wind.slice():[0,0,0]};}
 function moments(points){
  const n=points.length;if(!n)return {mean:null,covariance:null,rms:null};
  const mean=[points.reduce((s,p)=>s+p.y,0)/n,points.reduce((s,p)=>s+p.z,0)/n];
  if(n<2)return {mean,covariance:null,rms:null};
  const cov=[[0,0],[0,0]];for(const p of points){const dy=p.y-mean[0],dz=p.z-mean[1];cov[0][0]+=dy*dy;cov[0][1]+=dy*dz;cov[1][1]+=dz*dz;}
  cov[0][0]/=n-1;cov[0][1]/=n-1;cov[1][0]=cov[0][1];cov[1][1]/=n-1;
  return {mean,covariance:cov,rms:Math.sqrt(Math.max(0,cov[0][0]+cov[1][1]))};
 }
 function quantile(xs,p){const a=xs.slice().sort((a,b)=>a-b),x=p*(a.length-1),j=Math.floor(x);return a[j]+(a[Math.min(j+1,a.length-1)]-a[j])*(x-j);}
 function finiteSamplingInterval(points,seed,B){if(points.length<20)return null;const random=rng(seed),values=[];for(let j=0;j<B;j++)values.push(moments(Array.from({length:points.length},()=>points[Math.floor(random()*points.length)])).rms);return [quantile(values,.025),quantile(values,.975)];}
 function summaryFor(scheme,c,final){
  const samples=scheme.samples,N=samples.length,arrivals=samples.filter(s=>s.status==='crossed_goal_plane'),n=arrivals.length,points=arrivals.map(s=>({y:s.yGoal,z:s.zGoal})),mom=moments(points),statusCounts=Object.fromEntries(STATES.map(s=>[s,0]));
  for(const s of samples)statusCounts[s.status]=(statusCounts[s.status]||0)+1;
  const numeric=['invalid_input','fit_failed','integration_failed','outside_model_domain'].reduce((s,k)=>s+(statusCounts[k]||0),0);
  const nominal=scheme.nominal,rmse=nominal.fit?.rmse,normalizedRMSE=Number.isFinite(rmse)&&c.sigma>0?rmse/c.sigma:null;
  const fitAcceptable=!!nominal.fit?.ok&&Number.isFinite(rmse)&&(c.sigma>0?normalizedRMSE<=c.selection.maxNormalizedRMSE:rmse<=1e-6);
  const timeMean=n?arrivals.reduce((s,a)=>s+a.tGoal,0)/n:null,timeStd=n>1?Math.sqrt(arrivals.reduce((s,a)=>s+(a.tGoal-timeMean)**2,0)/(n-1)):null;
  const interval=final?finiteSamplingInterval(points,(c.seed+0x51d7348)>>>0,c.intervalReplicates):null;
  const exclusionReasons=[];
  if(c.sigma===0)exclusionReasons.push('σ=0：只有确定性结果，不产生统计密度或稳定性推荐');
  if(!fitAcceptable)exclusionReasons.push('名义观测残差未通过声明的容许门限');
  if(N<c.selection.minSamples||n<c.selection.minSamples)exclusionReasons.push('可用样本不足');
  if(N&&n/N<c.selection.minArrivalFraction)exclusionReasons.push('到达样本比例低，避免将截断后的窄分布误判为稳定');
  if(N&&numeric/N>c.selection.maxNumericalFailureFraction)exclusionReasons.push('拟合/输入/数值失败比例超过门限');
  if(!Number.isFinite(mom.rms))exclusionReasons.push('没有可比较的二维分散程度');
  const density=D.buildDensity(points,{...c.bandwidth,supportZMin:c.base.R,planeArrivalMass:N?{nArrival:n,nTotal:N}:null});
  if(c.sigma===0){density.status='degenerate';density.reason='zero_declared_input_uncertainty';}
  let kernelGoalMass=null;if(density.status==='ok')kernelGoalMass=D.integrateRectangle(density,{yMin:-c.base.goalHalfWidth,yMax:c.base.goalHalfWidth,zMin:0,zMax:c.base.goalHeight});
  const p=nominal.prediction,hasNominal=p?.status==='crossed_goal_plane';
  const marginal=n?{y:[quantile(arrivals.map(a=>a.yGoal),.025),quantile(arrivals.map(a=>a.yGoal),.975)],z:[quantile(arrivals.map(a=>a.zGoal),.025),quantile(arrivals.map(a=>a.zGoal),.975)],t:[quantile(arrivals.map(a=>a.tGoal),.025),quantile(arrivals.map(a=>a.tGoal),.975)],meaning:'各坐标各自2.5%–97.5%经验分位范围，不是同时95%置信区间或已校准真实覆盖'}:null;
  const summary={total:N,arrived:n,arrivalFraction:N?n/N:null,statusCounts,numericalFailureCount:numeric,numericalFailureFraction:N?numeric/N:null,boundaryHits:samples.filter(s=>s.boundaryHit).length,marginalRanges:marginal,
   normalizedRMSE,nominalRMSE_m:rmse??null,fitAcceptable,meanYZ:mom.mean,covarianceYZ:mom.covariance,stabilityVariance_m2:mom.rms===null?null:mom.rms**2,stabilityRMS_m:mom.rms,stabilityInterval_m:interval,
   intervalMeaning:'对到达样本的再抽样百分位区间；仅描述当前Monte Carlo分散度估计的有限样本变化，不是实球预测区间',arrivalTimeStd_s:timeStd,
   rmsDisplacementFromNominal_m:n&&hasNominal?Math.sqrt(arrivals.reduce((s,a)=>s+(a.yGoal-p.yGoal)**2+(a.zGoal-p.zGoal)**2,0)/n):null,
   eligible:exclusionReasons.length===0,exclusionReasons,empiricalCenterFractionOfArrivals:n?arrivals.filter(s=>s.centerInside).length/n:null,kernelGoalMass,
   probabilityMeaning:'基于声明观测扰动和拟合算法的条件输出；nArrival/N还受数值失败影响，不是实测到达概率'};
  return {summary,density};
 }
 function decide(schemes,c,final){
  const eligible=schemes.filter(s=>s.summary.eligible).sort((a,b)=>a.summary.stabilityRMS_m-b.summary.stabilityRMS_m);
  const policy={...c.selection,criterion:'sqrt(trace(Cov(y,z)))',kind:'predeclared operational stability rule, not proof of model truth',intervalRule:'Overlap of resampling intervals is a conservative ambiguity flag, not a significance test',boundaryRule:'Boundary counts are retained and warned; boundary-hit samples are not dropped'};
  if(!eligible.length)return {status:schemes.some(s=>s.summary.total<c.selection.minSamples)?'insufficient_samples':'no_eligible_scheme',preferredId:null,candidateIds:[],uniqueSupported:false,criterion:policy.criterion,policy,message:'尚无符合门限的可推荐方案；保留全部结果，不强行选胜者。'};
  const first=eligible[0],best=first.summary.stabilityInterval_m;
  const candidates=final&&best?eligible.filter(s=>{const ci=s.summary.stabilityInterval_m;return !ci||ci[0]<=best[1]&&ci[1]>=best[0];}):eligible;
  const unique=final&&candidates.length===1;
  return {status:unique?'conditional_minimum':'ambiguous',preferredId:first.id,candidateIds:candidates.map(s=>s.id),uniqueSupported:unique,criterion:policy.criterion,policy,
   message:unique?`${first.id}在当前噪声、门限和分散度指标下最稳定；这是条件性选择，不证明其受力机制为真实唯一模型。`:`当前最小数值来自${first.id}，但有限抽样区间不能清楚区分候选；暂不作唯一稳定推荐。`};
 }
 function compareSchemes(schemes){const pairs=[];for(let i=0;i<schemes.length;i++)for(let j=i+1;j<schemes.length;j++){
  const a=schemes[i],b=schemes[j],pa=a.nominal.prediction,pb=b.nominal.prediction,paired=a.samples.filter((s,k)=>s.status==='crossed_goal_plane'&&b.samples[k]?.status==='crossed_goal_plane');
  pairs.push({a:a.id,b:b.id,nominalSeparation_m:pa.status==='crossed_goal_plane'&&pb.status==='crossed_goal_plane'?Math.hypot(pa.yGoal-pb.yGoal,pa.zGoal-pb.zGoal):null,pairedArrivals:paired.length,
   pairedPredictionSeparationRMS_m:paired.length?Math.sqrt(paired.reduce((v,s)=>{const t=b.samples[s.index];return v+(s.yGoal-t.yGoal)**2+(s.zGoal-t.zGoal)**2;},0)/paired.length):null});}
  const nominal=schemes.filter(s=>s.nominal.prediction.status==='crossed_goal_plane').map(s=>({y:s.nominal.prediction.yGoal,z:s.nominal.prediction.zGoal}));
  return {pairs,betweenNominalSpread_m:moments(nominal).rms,interpretation:'方案间差异单独报告；候选组合没有自动概率权重，不将任意枚举方案的分散度当作现实总体概率。'};
 }
 function runStudy(input,onProgress){
  const began=Date.now(),config=normalize(input),issues=validate(config);if(issues.length)throw Error(issues.join('；'));
  const trials=makeTrials(config),schemes=config.schemes.map(s=>{const cfg=schemeConfig(config,s);return {...s,label:s.label||PRESETS[s.preset].label,physicsConfig:cfg,nominal:P.solve(cfg),samples:[]};});
  function snapshot(final){for(const s of schemes)Object.assign(s,summaryFor(s,config,final));return {schemes,decision:decide(schemes,config,final)};}
  if(onProgress){const snap=snapshot(false);onProgress({completed:0,total:trials.length,...snap,done:false});}
  for(let i=0;i<trials.length;i++){
   const trial=trials[i];
   for(const s of schemes){
    let run;try{run=P.solve({...s.physicsConfig,observations:trial.observations,fitStart:s.nominal.theta||undefined});}catch(error){run={prediction:{status:'integration_failed'},fit:{ok:false,boundaryHit:false,rmse:null},issues:[String(error)]};}
    const p=run.prediction;s.samples.push({index:i,trialSeed:trial.seed,status:p.status,yGoal:p.yGoal??null,zGoal:p.zGoal??null,tGoal:p.tGoal??null,eventTime:p.eventTime??null,centerInside:p.centerInside??false,clearanceInside:p.clearanceInside??false,theta:run.theta??null,constantAcceleration:run.fit?.constantAcceleration??null,fitOK:run.fit?.ok??false,rmse:run.fit?.rmse??null,boundaryHit:run.fit?.boundaryHit??false,issues:run.issues||[]});
   }
   if(onProgress&&((i+1)%10===0||i===trials.length-1)){const snap=snapshot(false);onProgress({completed:i+1,total:trials.length,schemes:snap.schemes.map(({samples,...s})=>s),decision:snap.decision,done:false});}
  }
  const snap=snapshot(true),assessmentSeed=(config.seed^0xa5b35705)>>>0,assessmentConfig={...config,seed:assessmentSeed,sampleCount:config.assessmentCount};
  const assessmentModels=config.assessmentEnabled?schemes.filter(s=>snap.decision.candidateIds.includes(s.id)).map(s=>({...s,samples:[]})):[],assessmentTrials=assessmentModels.length?makeTrials(assessmentConfig):[];
  for(let i=0;i<assessmentTrials.length;i++){
   const trial=assessmentTrials[i];
   for(const s of assessmentModels){let run;try{run=P.solve({...s.physicsConfig,observations:trial.observations,fitStart:s.nominal.theta||undefined});}catch(error){run={prediction:{status:'integration_failed'},fit:{ok:false},issues:[String(error)]};}
    const p=run.prediction;s.samples.push({index:i,trialSeed:trial.seed,status:p.status,yGoal:p.yGoal??null,zGoal:p.zGoal??null,tGoal:p.tGoal??null,theta:run.theta??null,constantAcceleration:run.fit?.constantAcceleration??null,fitOK:run.fit?.ok??false,rmse:run.fit?.rmse??null,boundaryHit:run.fit?.boundaryHit??false,centerInside:p.centerInside??false,issues:run.issues||[]});}
   if(onProgress&&((i+1)%10===0||i===assessmentTrials.length-1))onProgress({completed:config.sampleCount,total:config.sampleCount,stage:'assessment',assessmentCompleted:i+1,assessmentTotal:assessmentTrials.length,schemes:snap.schemes.map(({samples,...s})=>s),decision:snap.decision,done:false});
  }
  for(const s of assessmentModels)Object.assign(s,summaryFor(s,assessmentConfig,true));
  const assessment={enabled:config.assessmentEnabled,seed:assessmentSeed,trialCount:assessmentTrials.length,trials:assessmentTrials,models:assessmentModels,selectionFrozen:true,interpretation:'独立伪随机流重新扰动所选及并列候选；选择冻结，不用此批样本再次挑最小值。仅量化指定输入误差与模型下的不确定度，不含未知模型偏差。'};
  return {version:VERSION,config,trials,schemes:snap.schemes,decision:snap.decision,assessment,comparison:compareSchemes(schemes),elapsedMs:Date.now()-began,
   mappingDefinition:{input:'two observed positions, declared independent Gaussian coordinate sigma, an explicitly enumerated candidate set, model assumptions, finite sampling seed and budgets',randomization:'the SAME observations + epsilon are passed to every candidate per index; selection and assessment use separate pseudorandom streams',output:'prediction mapping F(p_0.05,p_0.10) -> candidate selection, nominal y/z/t predictions and independently assessed conditional uncertainty; KDE is a supplementary representation',function:'F(d; config) = Predict_mhat(Fit_mhat(d)), with mhat = argmin eligible_m sqrt(trace(Cov_pilot(y,z))); ties retained. f_m(y,z) = I[z>=R] sum_i w_i exp(-((y-y_i)^2+(z-z_i)^2)/(2h^2))/(2*pi*h^2*Phi((z_i-R)/h))',selectedFunction:snap.decision.preferredId,selectedFunctionIsProvisional:!snap.decision.uniqueSupported,modelWeightsInferred:false,physicalValidationCompleted:false,
   assumptions:['Initial position and timestamps exact','Independent, equal-variance Gaussian coordinate perturbations are user assumptions, not measured','Wind for wind preset is a specified scenario, not inferred from six coordinates','Effective Magnus b is not true angular velocity','Kernel bandwidth controls density smoothing and is not a fitted physical uncertainty','Each Gaussian kernel is truncated and normalized at ball-center ground height R; no truncation at goal posts or crossbar','Minimum variance is a conditional decision criterion; independent measurements are needed to establish predictive accuracy']}};
 }
 return {VERSION,defaultStudyConfig,generateSchemes,validate:input=>validate(normalize(input)),makeTrials:input=>makeTrials(normalize(input)),runStudy,schemeConfig,moments};
}
if(typeof module!=='undefined'&&module.exports)module.exports={createStudy};
if(typeof globalThis!=='undefined')globalThis.createStudy=createStudy;
