/* Independent implementation of rigid-sphere forward dynamics, M4.
 * Literature formula: Kiratidis & Leinweber, arXiv:1710.02784v2 (2018),
 * equations 5, 7, 9 and Table 1. No article, experimental image or source code
 * is redistributed here; the arXiv article carries its non-exclusive licence.
 * The linear lift law, oblique-axis extension, inertia and spin damping are
 * explicit modelling assumptions, not measured properties of the UPC ball.
 */
function createAdvanced(P) {
  'use strict';
  const VERSION='upc-m4-forward-1.0.0';
  const SOURCE='https://arxiv.org/html/1710.02784';
  const TABLE={
    Tango12:[.5452,12.86,1.304,.1657,.1953,16.22,35],
    Teamgeist:[.4927,12.58,1.071,.1440,.1540,23.17,35],
    Jabulani:[.4839,18.69,1.377,.1413,.1780,23.29,35],
    Brazuca:[.4740,12.92,1,.1657,.2112,14.61,35]
  };
  const norm=v=>Math.hypot(...v);
  const cross=(a,b)=>[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
  const finiteVector=(x,n)=>Array.isArray(x)&&x.length===n&&x.every(Number.isFinite);
  function defaultConfig(){return {...P.defaultConfig(),omega:[0,0,50],quaternion:[1,0,0,0],
    spinDecayTau:'infinite',kappa:2/3,dragModel:'constant',ballType:'Teamgeist',
    liftModel:'linear',liftSlope:.8,mu:1.81e-5,windShear:[0,0,0],gustAmplitude:[0,0,0],gustPeriod:1};}
  function config(input){return {...defaultConfig(),...(input||{})};}
  const advancedKeys=['omega','quaternion','spinDecayTau','kappa','dragModel','ballType','liftModel','liftSlope','mu','windShear','gustAmplitude','gustPeriod'];
  function baseInput(input){const c={...(input||{})};for(const k of advancedKeys)delete c[k];return c;}
  const tau=c=>c.spinDecayTau==='infinite'||c.spinDecayTau===Infinity?Infinity:c.spinDecayTau;
  function validateConfig(input){
    const c=config(input);if(c.model!=='M4')return P.validateConfig(baseInput(c));
    // P0 validates shared SI quantities, times, geometry and step bounds.
    const issues=P.validateConfig({...baseInput(c),model:'M3'});
    if(c.mode!=='forward')issues.push('M4只支持给定真实角速度的正向研究，未实现两帧自旋反演。');
    if(!finiteVector(c.windShear,3)||norm(c.windShear)>100)issues.push('风切变须为有限三分量，数值工作范围≤100s⁻¹。');
    if(!finiteVector(c.gustAmplitude,3)||norm(c.gustAmplitude)>100)issues.push('阵风幅度须为有限三分量，范围≤100m/s。');
    if(!Number.isFinite(c.gustPeriod)||c.gustPeriod<.01)issues.push('确定性阵风周期须至少0.01s（数值工作范围）。');
    if(!finiteVector(c.omega,3)||norm(c.omega)>2000)issues.push('omega须为3个有限rad/s数值，数值工作范围不超过2000 rad/s。');
    if(!finiteVector(c.quaternion,4)||Math.abs(norm(c.quaternion)-1)>1e-8)issues.push('quaternion须为单位四元数[w,x,y,z]，不得静默修正输入。');
    if(!Number.isFinite(c.kappa)||c.kappa<=0||c.kappa>1)issues.push('惯量系数kappa须在(0,1]，默认球壳2/3是明确假设。');
    if(!(tau(c)===Infinity||(Number.isFinite(tau(c))&&tau(c)>=1e-4)))issues.push('spinDecayTau须为infinite或至少0.0001s的有限正数（数值工作范围）。');
    if(!Number.isFinite(c.mu)||c.mu<=0)issues.push('动力黏度mu必须为正有限Pa·s。');
    if(!Number.isFinite(c.liftSlope)||c.liftSlope<0||c.liftSlope>10)issues.push('线性升力斜率须在0–10的声明情景范围。');
    if(!['constant','paper2018'].includes(c.dragModel))issues.push('未知M4阻力关系。');
    if(!Object.prototype.hasOwnProperty.call(TABLE,c.ballType))issues.push('球型须为Tango12/Teamgeist/Jabulani/Brazuca。');
    if(!['linear','teamgeist2018'].includes(c.liftModel))issues.push('未知M4升力关系。');
    if(c.liftModel==='teamgeist2018'&&(c.ballType!=='Teamgeist'||c.dragModel!=='paper2018'))issues.push('Teamgeist2018升力须同时选择Teamgeist与paper2018阻力。');
    return issues;
  }
  function baseCd(U,ballType){
    const [a,vc,vs,bmin,bmax,vmin,vmax]=TABLE[ballType];
    return (a-bmin)/(1+Math.exp((U-vc)/vs))+bmin+
      (U-vmin)/(1+Math.exp(-(U-vmin)/vs))*(bmax-bmin)/(vmax-vmin);
  }
  function coefficients(U,omega,w,c){
    const om=norm(omega),wx=U>0?w.map(x=>x/U):[0,0,0];
    const omegaPerp=norm(cross(omega,wx));
    const S=U>0?c.R*om/U:null, Sperp=U>0?c.R*omegaPerp/U:null;
    const sinAngle=om>0?omegaPerp/om:0;
    let Cd=c.Cd,CL=0,status='ok',message='';
    if(c.dragModel==='paper2018'){
      Cd=baseCd(U,c.ballType);
      if(U<7||U>35){status='outside_model_domain';message='paper2018阻力超出项目选定7–35m/s保护域；该区间不是原文置信区间。';}
    }
    if(c.liftModel==='teamgeist2018'){
      // Experimental spin drag slope, plus the paper's proposed lift interpolation.
      // Applying these transverse-spin relations at oblique angles is an assumption.
      if(S!==null)Cd+=.514*S;
      if(U<23.17||U>35||S===null||S>.4){status='outside_model_domain';message='Teamgeist情景要求23.17≤U≤35m/s且S≤0.4；不作域外外推。';}
      if(S!==null)CL=1.15*Math.pow(S,.83)*(.4927-Math.min(Cd,.4927))/(.4927-(.144+.514*.19));
    } else if(Sperp!==null){
      CL=c.liftSlope*Sperp;
    }
    return {Cd,CL,S,S_perp:Sperp,sinAngle,Re:c.rho*U*2*c.R/c.mu,status,message};
  }
  function forceM4(t,state,input){
    const c=config(input),s=Array.isArray(state)?state:[...state.position,...state.velocity,...(state.omega||c.omega),...(state.q||c.quaternion)];
    const omega=s.length>=9?s.slice(6,9):c.omega;
    const airVelocity=c.wind.map((x,i)=>x+c.windShear[i]*s[2]+c.gustAmplitude[i]*Math.sin(2*Math.PI*t/c.gustPeriod));
    const w=s.slice(3,6).map((x,i)=>x-airVelocity[i]),U=norm(w);
    const coeff=coefficients(U,omega,w,c),A=Math.PI*c.R*c.R;
    const drag=w.map(x=>-.5*c.rho*A*coeff.Cd*U*x);
    const direction=cross(omega,w),dn=norm(direction),scale=.5*c.rho*A*U*U*coeff.CL;
    const magnus=direction.map(x=>dn>0?scale*x/dn*(c.liftModel==='teamgeist2018'?coeff.sinAngle:1):0);
    const gravity=[0,0,-c.m*c.g],total=gravity.map((x,i)=>x+drag[i]+magnus[i]);
    const inertia=c.kappa*c.m*c.R*c.R,decay=tau(c);
    const torque=omega.map(x=>decay===Infinity?0:-inertia*x/decay);
    return {gravity,drag,magnus,total,acceleration:total.map(x=>x/c.m),relativeAirVelocity:w,airVelocity,airSpeed:U,
      effectiveB:null,Cd:coeff.Cd,CL:coeff.CL,CL_effective:U>0?norm(magnus)/(.5*c.rho*A*U*U):null,
      Re:coeff.Re,S:coeff.S,S_perp:coeff.S_perp,omega:omega.slice(),torque,inertia,
      domainStatus:coeff.status,domainMessage:coeff.message,
      forceInterpretation:'M4真实角速度正向情景；升力角度推广、惯量及自旋衰减均须按假设解释。'};
  }
  function evaluateForces(t,state,input,theta){
    return input&&input.model==='M4'?forceM4(t,state,input):P.evaluateForces(t,state,baseInput(input),theta);
  }
  function derivative(t,s,c){
    const f=forceM4(t,s,c),o=s.slice(6,9),q=s.slice(9,13),T=tau(c);
    const oq=[-o[0]*q[1]-o[1]*q[2]-o[2]*q[3],
      o[0]*q[0]+o[1]*q[3]-o[2]*q[2],
      -o[0]*q[3]+o[1]*q[0]+o[2]*q[1],
      o[0]*q[2]-o[1]*q[1]+o[2]*q[0]];
    return [...s.slice(3,6),...f.acceleration,...o.map(x=>T===Infinity?0:-x/T),...oq.map(x=>x/2)];
  }
  function step(t,s,h,c){
    const k1=derivative(t,s,c),k2=derivative(t+h/2,s.map((x,i)=>x+h*k1[i]/2),c),
      k3=derivative(t+h/2,s.map((x,i)=>x+h*k2[i]/2),c),k4=derivative(t+h,s.map((x,i)=>x+h*k3[i]),c);
    const next=s.map((x,i)=>x+h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6);
    const qn=norm(next.slice(9,13));
    if(qn>0)for(let i=9;i<13;i++)next[i]/=qn;
    return next;
  }
  function row(t,s,c){const forces=forceM4(t,s,c);return {t,position:s.slice(0,3),velocity:s.slice(3,6),
    acceleration:forces.acceleration,forces,omega:s.slice(6,9),q:s.slice(9,13),Cd:forces.Cd,CL:forces.CL,
    Re:forces.Re,S:forces.S,torque:forces.torque};}
  function blank(status,message=''){return {status,yGoal:null,zGoal:null,tGoal:null,eventTime:null,
    centerInside:false,clearanceInside:false,message};}
  function integrateTrajectory(input,theta,options){
    if(!input||input.model!=='M4')return P.integrateTrajectory(baseInput(input),theta,options);
    const c=config(input),issues=validateConfig(c),collect=!options||options.collect!==false;
    if(issues.length)return {trajectory:[],prediction:blank('invalid_input',issues.join('; ')),issues};
    const v=theta?theta.slice(0,3):c.velocity;
    if(!finiteVector(v,3)||norm(v)>100)return {trajectory:[],prediction:blank('invalid_input','实际初速超出100m/s数值工作范围。')};
    let s=[...c.initial,...v,...c.omega,...c.quaternion],t=0;
    const trajectory=collect?[row(0,s,c)]:[],initialForce=forceM4(0,s,c);
    let pred=blank('no_crossing_within_horizon');
    if(s[2]<=c.R+1e-12&&(s[5]<0||(s[5]===0&&initialForce.acceleration[2]<=0))){
      pred=blank('ground_before_goal');pred.eventTime=0;return {trajectory,prediction:pred};
    }
    if(initialForce.domainStatus!=='ok'){
      pred=blank('outside_model_domain',initialForce.domainMessage);pred.eventTime=0;return {trajectory,prediction:pred};
    }
    let count=0;
    for(;count<200000&&t<c.horizon-1e-13;count++){
      const dt=Math.min(c.h,c.horizon-t,tau(c)/20,c.gustPeriod/40,.1/Math.max(norm(s.slice(6,9)),1));
      const next=step(t,s,dt,c);
      if(next.some(x=>!Number.isFinite(x)||Math.abs(x)>1e9)){
        pred=blank('integration_failed','非有限状态或数值范围溢出');pred.eventTime=t;break;
      }
      const events=[];
      if(s[0]<c.goalX&&next[0]>=c.goalX)events.push({kind:'goal',index:0,target:c.goalX});
      if(next[2]<=c.R&&(s[2]>c.R||(t===0&&s[5]>0)))events.push({kind:'ground',index:2,target:c.R});
      const nextForce=forceM4(t+dt,next,c);
      if(nextForce.domainStatus!=='ok')events.push({kind:'domain'});
      for(const e of events){
        let lo=0,hi=dt;
        for(let j=0;j<45&&hi-lo>1e-12;j++){
          const mid=(lo+hi)/2,q=step(t,s,mid,c);
          const after=e.kind==='domain'?forceM4(t+mid,q,c).domainStatus!=='ok':e.kind==='goal'?q[0]>=e.target:q[2]<=e.target;
          if(after)hi=mid;else lo=mid;
        }
        e.dt=(lo+hi)/2;e.state=step(t,s,e.dt,c);
      }
      if(events.length){
        events.sort((a,b)=>a.dt-b.dt);const e=events[0];t+=e.dt;s=e.state;
        if(e.kind==='domain')pred=blank('outside_model_domain',nextForce.domainMessage);
        else if(e.kind==='ground'){s[2]=c.R;pred=blank('ground_before_goal');}
        else{
          s[0]=c.goalX;pred=blank('crossed_goal_plane');pred.yGoal=s[1];pred.zGoal=s[2];pred.tGoal=t;
          pred.centerInside=Math.abs(s[1])<c.goalHalfWidth&&s[2]>0&&s[2]<c.goalHeight;
          pred.clearanceInside=Math.abs(s[1])<c.goalHalfWidth-c.R&&s[2]>=c.R&&s[2]<c.goalHeight-c.R;
        }
        pred.eventTime=t;if(collect)trajectory.push(row(t,s,c));break;
      }
      t+=dt;s=next;if(collect)trajectory.push(row(t,s,c));
    }
    if(pred.eventTime===null){pred.eventTime=t;if(count>=200000)pred= {...blank('integration_failed','超过200000内部步预算'),eventTime:t};}
    return {trajectory,prediction:pred,diagnostics:{internalSteps:count,quaternionConvention:'scalar-first, body-to-world, qdot=0.5*[0,omega_world]⊗q',
      numericalMethod:'RK4; step=min(h,tau/20,gustPeriod/40,0.1/max(|omega|,1)); event bisection; unit-quaternion projection'}};
  }
  function predictPositions(input,theta,times){
    if(!input||input.model!=='M4')return P.predictPositions(baseInput(input),theta,times);
    const c=config(input),issues=validateConfig(c);
    if(issues.length||!Array.isArray(times)||times.some((x,i)=>!Number.isFinite(x)||x<0||(i&&x<times[i-1])))throw new Error('M4位置采样输入无效');
    const run=integrateTrajectory(c,theta),rows=run.trajectory;
    return times.map(target=>{
      if(!rows.length||target>run.prediction.eventTime+1e-10)return null;
      let lo=0,hi=rows.length-1;while(lo<hi){const mid=Math.ceil((lo+hi)/2);if(rows[mid].t<=target)lo=mid;else hi=mid-1;}
      const r=rows[lo],s=[...r.position,...r.velocity,...r.omega,...r.q];
      return step(r.t,s,Math.max(0,target-r.t),c).slice(0,3);
    });
  }
  function solve(input){
    if(!input||input.model!=='M4'){const result=P.solve(baseInput(input));return {...result,config:{...config(input),...result.config}};}
    const started=Date.now(),c=config(input),issues=validateConfig(c),theta=finiteVector(c.velocity,3)?[...c.velocity,c.Cd,0,0]:null;
    const fit={ok:issues.length===0,rmse:null,residuals:[],starts:[],boundaryHit:false,
      diagnostics:{algorithm:'forward only',freeParameters:[],message:'真实omega由用户给定，未由两帧反演；有效b在M4不参与方程。'}};
    const result=integrateTrajectory(c,theta);
    return {version:VERSION,config:c,theta,fit,issues,...result,
      predictedObservations:issues.length?[]:predictPositions(c,theta,c.observations.map(o=>o.t)),elapsedMs:Date.now()-started,
      assumptions:['真实角速度初值是用户情景，非两帧识别值；b不参与M4',
        '各向同性球壳惯量I=kappa*m*R²；默认kappa=2/3是假设',
        '空气力矩=-I*omega/tau，自旋指数衰减是假设；infinite表示无力矩',
        c.liftModel==='linear'?'CL=liftSlope*R*|omega_perp|/U，默认斜率0.8是可调假设':'Teamgeist2018采用实验拟合及论文物理插值；斜轴sin角度推广为新增假设',
        c.dragModel==='paper2018'?'非旋转Cd拟合仅对应所选球型；7–35m/s为项目保护域，不是实验置信区间':'常Cd是指定情景值，未标定真实球型',
        '风=均匀风+高度线性切变+正弦阵风；属于单向给定外界风，非CFD',
        '未完成独立实球实验确认；未求解双向CFD或接缝姿态依赖'],
      equation:'r′=v; m v′=m g+F_D+F_L; I omega′=tau_air; q′=0.5[0,omega_world]⊗q',
      sources:[{title:'Kiratidis & Leinweber (2018), Eqs.5/7/9 and Table1',url:SOURCE,
        licence:'arXiv non-exclusive article licence; article/code/images not redistributed; independently implemented equations'}]};
  }
  function fitObservations(input){
    if(input&&input.model==='M4')return {ok:false,theta:null,rmse:null,residuals:[],starts:[],boundaryHit:false,
      diagnostics:{issues:['M4尚未实现反演；不能仅靠两帧唯一恢复真实自旋。']}};
    return P.fitObservations(baseInput(input));
  }
  function uncertainty(input,...args){if(input&&input.model==='M4')throw new Error('M4未实现概率推断；请使用正向参数情景。');return P.uncertainty(baseInput(input),...args);}
  function computeSample(input,...args){if(input&&input.model==='M4')throw new Error('M4不支持两帧bootstrap。');return P.computeSample(baseInput(input),...args);}
  function compare(input,comparison){if(input&&input.model==='M4')throw new Error('M4暂以保存配置比较正向情景，不自动转换真实omega为有效b。');return P.compare(baseInput(input),comparison);}
  return {...P,ADVANCED_VERSION:VERSION,defaultConfig,defaultAdvancedConfig:defaultConfig,validateConfig,
    solve,run:solve,evaluateForces,integrateTrajectory,predictPositions,fitObservations,uncertainty,computeSample,compare,
    advancedDragTable:()=>JSON.parse(JSON.stringify(TABLE)),advancedSource:SOURCE};
}
if(typeof globalThis!=='undefined'){
  globalThis.createAdvanced=createAdvanced;
  if(globalThis.UPC)globalThis.UPC=createAdvanced(globalThis.UPC);
}
if(typeof module!=='undefined'&&module.exports){module.exports=createAdvanced;module.exports.createAdvanced=createAdvanced;}
