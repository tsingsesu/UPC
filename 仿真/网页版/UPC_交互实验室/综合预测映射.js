'use strict';
// SI units. F([x05,y05,z05],[x10,y10,z10]) reruns inverse fits and model selection.
const createPhysics=function createPhysics() {
  'use strict';
  const VERSION = 'upc-field-combination-js-1.1.0';
  const LB = [10, -20, -5, .05, -1.5, -1.5];
  const UB = [50, 20, 25, .60, 1.5, 1.5];
  const ADVANCED_FIELDS = ['omega','quaternion','spinDecayTau','kappa','dragModel','ballType','liftModel','liftSlope','mu','windShear','gustAmplitude','gustPeriod'];
  const clone = x => JSON.parse(JSON.stringify(x));
  const dot = (a,b) => a.reduce((s,x,i)=>s+x*b[i],0);
  const norm = a => Math.hypot(...a);
  const cross = (a,b) => [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
  const finite3 = x => Array.isArray(x)&&x.length===3&&x.every(Number.isFinite);
  const clamp = (x,l,u)=>Math.min(u,Math.max(l,x));
  function defaultConfig() {
    return {schemaVersion:'1.0',model:'M3',mode:'inverse',initial:[0,0,.11],
      observations:[{t:.05,position:[1.5,.45,.44]},{t:.1,position:[2.96,.86,.74]}],
      m:.43,R:.11,rho:1.225,g:9.80665,wind:[0,0,0],sigma:.01,
      h:.001,fitH:.005,horizon:2,goalX:11,goalHalfWidth:3.66,goalHeight:2.44,
      velocity:[30.4,9.4,6.9],Cd:.3,b:[0,0,0],constantAcceleration:[-16,-16,-12],forceMask:{drag:true,magnus:true}};
  }
  function normalize(input) { return Object.assign(defaultConfig(),clone(input||{})); }
  function validateConfig(input) {
    const c=normalize(input), issues=[];
    const allowed=[...Object.keys(defaultConfig()),...ADVANCED_FIELDS,'fitStart','airPropertiesProvenance','metadata','extensions','scenarioId','provenance'];
    for(const key of Object.keys(c))if(!allowed.includes(key))issues.push('未知配置字段 '+key+'；元信息请放metadata/extensions，不参与计算');
    if(c.schemaVersion!=='1.0') issues.push('不支持的配置版本');
    if(!['M0','M1','M2','M3'].includes(c.model)) issues.push('模型必须为 M0–M3');
    if(!['forward','inverse'].includes(c.mode)) issues.push('模式必须为 forward/inverse');
    if(!c.forceMask||typeof c.forceMask!=='object'||Array.isArray(c.forceMask)||
      typeof c.forceMask.drag!=='boolean'||typeof c.forceMask.magnus!=='boolean'||
      Object.keys(c.forceMask).some(k=>!['drag','magnus'].includes(k)))issues.push('forceMask须仅含drag和magnus两个布尔开关；M0–M2中不参与方程');
    for(const k of ['initial','wind','velocity','b','constantAcceleration']) if(!finite3(c[k])) issues.push(k+' 必须是三个有限数值');
    for(const k of ['m','R','rho','h','fitH','horizon','goalHalfWidth','goalHeight']) if(!Number.isFinite(c[k])||c[k]<=0) issues.push(k+' 必须为正有限数');
    for(const k of ['g','sigma','Cd']) if(!Number.isFinite(c[k])||c[k]<0) issues.push(k+' 必须是非负有限数');
    if(!Number.isFinite(c.goalX)) issues.push('goalX 无效');
    if(c.h<.00001||c.h>.02||c.fitH<.00001||c.fitH>.02) issues.push('积分步长须在 0.00001–0.02 s 内');
    if(c.horizon>20) issues.push('首飞时间上限不得超过20 s');
    if(finite3(c.initial)&&(c.initial[2]<c.R-1e-12||c.initial[0]>=c.goalX)) issues.push('初始球心低于R或不在球门平面前方');
    if(!Array.isArray(c.observations)||c.observations.length!==2) issues.push('必须有两帧观测');
    else { let prev=0; for(const o of c.observations) {
      if(!o||!Number.isFinite(o.t)||o.t<=prev||o.t>c.horizon||!finite3(o.position)) issues.push('观测时刻须递增且坐标有限');
      else if(c.mode==='inverse'&&o.position[2]<c.R-1e-12) issues.push('观测球心低于接触高度R');
      if(o) prev=o.t;
    }}
    if(finite3(c.velocity)&&norm(c.velocity)>100) issues.push('速度超出本程序100 m/s数值工作范围');
    if(finite3(c.wind)&&norm(c.wind)>100) issues.push('风速超出本程序100 m/s数值工作范围');
    if(finite3(c.b)&&norm(c.b)>10) issues.push('有效b超出本程序10 s⁻¹工作范围');
    if(c.Cd>2) issues.push('Cd超出本程序0–2工作范围');
    if(finite3(c.constantAcceleration)&&norm(c.constantAcceleration)>500) issues.push('常加速度超出500 m/s²工作范围');
    return issues;
  }
  function forceFlags(c){return {drag:c.model==='M2'||(c.model==='M3'&&c.forceMask?.drag===true),magnus:c.model==='M3'&&c.forceMask?.magnus===true};}
  function activeIndices(c){const f=forceFlags(c);return [0,1,2,...(f.drag?[3]:[]),...(f.magnus?[4,5]:[])];}
  function effectiveTheta(c,th){const out=th.slice();if(out.length!==6)return out;const f=forceFlags(c);if(!f.drag)out[3]=0;if(!f.magnus){out[4]=0;out[5]=0;}return out;}
  function thetaFromConfig(c) {return effectiveTheta(c,[...c.velocity,c.Cd,c.b[1],c.b[2]]);}
  function acceleration(s,c,th) {
    if(c.model==='M1') return c.constantAcceleration.slice();
    if(c.model==='M0') return [0,0,-c.g];
    const w=[s[3]-c.wind[0],s[4]-c.wind[1],s[5]-c.wind[2]];
    const f=forceFlags(c),d=f.drag?-c.rho*Math.PI*c.R*c.R/(2*c.m)*th[3]*norm(w):0;
    const b=f.magnus?[c.mode==='inverse'?0:c.b[0],th[4],th[5]]:[0,0,0];
    const lift=cross(b,w); return [d*w[0]+lift[0],d*w[1]+lift[1],-c.g+d*w[2]+lift[2]];
  }
  function evaluateForces(t,state,input,theta) {
    const c=normalize(input),s=Array.isArray(state)?state:[...state.position,...state.velocity],th=effectiveTheta(c,theta||thetaFromConfig(c));
    const w=[s[3]-c.wind[0],s[4]-c.wind[1],s[5]-c.wind[2]],U=norm(w);
    const flags=forceFlags(c),useDrag=flags.drag;
    const drag=w.map(x=>useDrag?-.5*c.rho*Math.PI*c.R*c.R*th[3]*U*x:0);
    const b=flags.magnus?[c.mode==='inverse'?0:c.b[0],th[4],th[5]]:[0,0,0];
    const magnus=cross(b,w).map(x=>c.m*x),gravity=[0,0,-c.m*c.g];
    const a=acceleration(s,c,th),total=a.map(x=>c.m*x);
    return {gravity,drag,magnus,total,acceleration:a,relativeAirVelocity:w,airSpeed:U,
      effectiveB:b,Cd:useDrag?th[3]:0,CL_effective:U>1e-12?2*norm(magnus)/(c.rho*Math.PI*c.R*c.R*U*U):null,
      forceInterpretation:c.model==='M1'?'常加速度代数延拓；total是m*a，非已分解气动力':'重力、平方阻力与有效Magnus力，单位N'};
  }
  function deriv(s,c,th) {return [s[3],s[4],s[5],...acceleration(s,c,th)];}
  function step(s,dt,c,th) {
    if(dt===0) return s.slice();
    const k1=deriv(s,c,th),k2=deriv(s.map((x,i)=>x+dt*k1[i]/2),c,th),
      k3=deriv(s.map((x,i)=>x+dt*k2[i]/2),c,th),k4=deriv(s.map((x,i)=>x+dt*k3[i]),c,th);
    return s.map((x,i)=>x+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6);
  }
  function predictPositions(input,theta,times) {
    const c=normalize(input);
    if(!Array.isArray(theta)||theta.length!==6||!theta.every(Number.isFinite)||!Array.isArray(times)||
      times.some((x,i)=>!Number.isFinite(x)||x<0||x>20||(i>0&&x<times[i-1]))||
      !Number.isFinite(c.fitH)||c.fitH<.00001||c.fitH>.02)throw new Error('观测预测参数或时间网格无效');
    theta=effectiveTheta(c,theta);let s=[...c.initial,...theta.slice(0,3)],t=0;
    return times.map(target=>{while(t<target-1e-13){const dt=Math.min(c.fitH,target-t);s=step(s,dt,c,theta);t+=dt;}return s.slice(0,3);});
  }
  function row(t,s,c,th) {const forces=evaluateForces(t,s,c,th);return {t,position:s.slice(0,3),velocity:s.slice(3),acceleration:forces.acceleration,forces};}
  function blankPrediction(status) {return {status,yGoal:null,zGoal:null,tGoal:null,eventTime:null,centerInside:false,clearanceInside:false};}
  function integrateTrajectory(input,theta,options) {
    const c=normalize(input),rawTheta=theta||thetaFromConfig(c),collect=!options||options.collect!==false;
    const th=Array.isArray(rawTheta)?effectiveTheta(c,rawTheta):rawTheta;
    const issues=validateConfig(c);if(!Array.isArray(rawTheta)||rawTheta.length!==6||!rawTheta.every(Number.isFinite)||rawTheta[3]<0)issues.push('theta须为六个有限数值且Cd非负');
    if(issues.length) return {trajectory:[],prediction:blankPrediction('invalid_input'),issues};
    if(norm(th.slice(0,3))>100) return {trajectory:[],prediction:blankPrediction('outside_model_domain'),issues:['实际初速超出100 m/s数值工作范围']};
    let s=[...c.initial,...th.slice(0,3)],t=0,pred=blankPrediction('no_crossing_within_horizon');
    const trajectory=collect?[row(0,s,c,th)]:[];
    const a0=acceleration(s,c,th);
    if(s[2]<=c.R+1e-12&&(s[5]<0||(s[5]===0&&a0[2]<=0))) {
      pred.status='ground_before_goal';pred.eventTime=0;return {trajectory,prediction:pred};
    }
    const maxSteps=Math.ceil(c.horizon/c.h)+2;
    for(let count=0;count<maxSteps&&t<c.horizon-1e-13;count++) {
      const dt=Math.min(c.h,c.horizon-t),next=step(s,dt,c,th);
      if(next.some(x=>!Number.isFinite(x)||Math.abs(x)>1e9)) {pred.status='integration_failed';pred.eventTime=t;break;}
      const events=[];
      if(s[0]<c.goalX&&next[0]>=c.goalX) events.push({i:0,target:c.goalX,status:'crossed_goal_plane'});
      if(next[2]<=c.R&&((s[2]>c.R)||(t===0&&s[5]>0))) events.push({i:2,target:c.R,status:'ground_before_goal'});
      for(const e of events) {
        let lo=0,hi=dt;
        for(let j=0;j<45&&hi-lo>1e-12;j++) {const mid=(lo+hi)/2,q=step(s,mid,c,th);
          if(e.i===0?q[e.i]>=e.target:q[e.i]<=e.target)hi=mid;else lo=mid;}
        e.dt=(lo+hi)/2;e.state=step(s,e.dt,c,th);
      }
      if(events.length){events.sort((a,b)=>a.dt-b.dt);const e=events[0];t+=e.dt;s=e.state;s[e.i]=e.target;
        pred.status=e.status;pred.eventTime=t;
        if(e.i===0){pred.yGoal=s[1];pred.zGoal=s[2];pred.tGoal=t;pred.normalVelocity=s[3];pred.nearTangent=Math.abs(s[3])<1e-3;
          if(pred.nearTangent)pred.warning='几乎与球门平面相切；到达时间对扰动敏感，不作高精度物理解读';
          pred.centerInside=Math.abs(s[1])<c.goalHalfWidth&&s[2]>0&&s[2]<c.goalHeight;
          pred.clearanceInside=Math.abs(s[1])<c.goalHalfWidth-c.R&&s[2]>=c.R&&s[2]<c.goalHeight-c.R;}
        if(collect)trajectory.push(row(t,s,c,th));break;
      }
      t+=dt;s=next;if(collect)trajectory.push(row(t,s,c,th));
    }
    if(pred.eventTime===null)pred.eventTime=t;
    return {trajectory,prediction:pred};
  }
  function linearSolve(A,b) {
    const n=b.length,M=A.map((r,i)=>[...r,b[i]]);
    for(let k=0;k<n;k++){let p=k;for(let i=k+1;i<n;i++)if(Math.abs(M[i][k])>Math.abs(M[p][k]))p=i;
      if(Math.abs(M[p][k])<1e-24)return null;[M[k],M[p]]=[M[p],M[k]];
      const d=M[k][k];for(let j=k;j<=n;j++)M[k][j]/=d;
      for(let i=0;i<n;i++)if(i!==k){const f=M[i][k];for(let j=k;j<=n;j++)M[i][j]-=f*M[k][j];}}
    return M.map(r=>r[n]);
  }
  function eigenSym(A) {
    const n=A.length,a=A.map(r=>r.slice()),v=Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
    for(let sweep=0;sweep<80*n*n;sweep++){let p=0,q=1,big=0;
      for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(a[i][j])>big){big=Math.abs(a[i][j]);p=i;q=j;}
      if(big<1e-16)break;const phi=.5*Math.atan2(2*a[p][q],a[q][q]-a[p][p]),c=Math.cos(phi),s=Math.sin(phi);
      const app=a[p][p],aqq=a[q][q],apq=a[p][q];a[p][p]=c*c*app-2*s*c*apq+s*s*aqq;a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;a[p][q]=a[q][p]=0;
      for(let k=0;k<n;k++)if(k!==p&&k!==q){const akp=a[k][p],akq=a[k][q];a[k][p]=a[p][k]=c*akp-s*akq;a[k][q]=a[q][k]=s*akp+c*akq;}
      for(let k=0;k<n;k++){const vp=v[k][p],vq=v[k][q];v[k][p]=c*vp-s*vq;v[k][q]=s*vp+c*vq;}}
    return a.map((r,i)=>({value:r[i],vector:v.map(row=>row[i])})).sort((a,b)=>b.value-a.value);
  }
  function jacobian(fun,q) {
    const base=fun(q),J=base.map(()=>Array(q.length).fill(0));
    for(let j=0;j<q.length;j++){const d=1e-5,qp=q.slice(),qm=q.slice();qp[j]=Math.min(1,q[j]+d);qm[j]=Math.max(0,q[j]-d);
      const p=fun(qp),m=fun(qm),den=qp[j]-qm[j];for(let i=0;i<base.length;i++)J[i][j]=(p[i]-m[i])/den;}
    return J;
  }
  function normal(J,r) {const n=J[0].length,A=Array.from({length:n},()=>Array(n).fill(0)),g=Array(n).fill(0);
    for(let i=0;i<J.length;i++)for(let j=0;j<n;j++){g[j]+=J[i][j]*r[i];for(let k=0;k<n;k++)A[j][k]+=J[i][j]*J[i][k];}return {A,g};}
  function fitObservationsRaw(input) {
    const c=normalize(input),issues=validateConfig(c);
    if(issues.length)return {ok:false,theta:null,rmse:null,residuals:[],starts:[],boundaryHit:false,diagnostics:{issues}};
    const times=c.observations.map(o=>o.t),obs=c.observations.map(o=>o.position),target=obs.flat();
    const indices=activeIndices(c),n=indices.length,parameterNames=['vx','vy','vz','Cd','by','bz'];
    function residual(th){return predictPositions(c,th,times).flat().map((x,i)=>x-target[i]);}
    if(c.model==='M0'||c.model==='M1'){
      let v,a=[0,0,-c.g];
      if(c.model==='M0'){const denom=dot(times,times);v=[0,1,2].map(j=>times.reduce((s,t,i)=>s+t*(obs[i][j]-c.initial[j]-.5*a[j]*t*t),0)/denom);}
      else {const t1=times[0],t2=times[1];a=[0,1,2].map(j=>2*((obs[1][j]-c.initial[j])/t2-(obs[0][j]-c.initial[j])/t1)/(t2-t1));v=[0,1,2].map(j=>(obs[0][j]-c.initial[j])/t1-.5*a[j]*t1);c.constantAcceleration=a;}
      const th=[...v,0,0,0],r=predictPositions(c,th,times).flat().map((x,i)=>x-target[i]);
      const eventResult=integrateTrajectory(c,th,{collect:false}),event=eventResult.prediction;
      const beforeObservation=event.status==='ground_before_goal'&&event.eventTime<times[1]-1e-9;
      const ok=!['invalid_input','integration_failed','outside_model_domain'].includes(event.status)&&!beforeObservation;
      const rejection=ok?null:beforeObservation?'ground_before_observation':event.status;
      return {ok,theta:th,constantAcceleration:c.model==='M1'?a:null,rmse:Math.sqrt(dot(r,r)/6),sse:dot(r,r),residuals:r,
        starts:[{theta:th,sse:dot(r,r),converged:true,accepted:ok,iterations:0}],boundaryHit:false,
        diagnostics:{algorithm:c.model==='M0'?'analytic linear least squares':'analytic constant-acceleration interpolation',rank:c.model==='M0'?3:6,condition:null,freeParameters:c.model==='M0'?['vx','vy','vz']:['vx','vy','vz','ax','ay','az'],rejectionReason:rejection,eventStatus:event.status,issues:eventResult.issues||[],message:ok?'解析基线；不证明物理真实性':beforeObservation?'拟合轨迹在观测前触地':'解析拟合未通过首飞/数值工作域检查：'+event.status+(eventResult.issues?.length?'；'+eventResult.issues.join('；'):'')}};
    }
    const scales=UB.map((u,i)=>u-LB[i]);
    const unpack=q=>{const th=Array(6).fill(0);indices.forEach((j,i)=>{th[j]=LB[j]+scales[j]*q[i];});return th;};
    const fun=q=>residual(unpack(q));
    const starts=[[24,3,5,.25,0,0],[30,6,8,.3,0,-.3],[20,0,3,.2,.2,.2],[40,10,12,.4,-.2,.4]];
    if(Array.isArray(c.fitStart)&&c.fitStart.length===6&&c.fitStart.every(Number.isFinite))starts.unshift(c.fitStart);
    const results=[];
    for(const st of starts){let q=indices.map(j=>clamp((st[j]-LB[j])/scales[j],0,1));let r=fun(q),cost=dot(r,r),lambda=1e-5,converged=false,reason='iteration_limit',iterations=0;
      for(iterations=0;iterations<180;iterations++){
        const J=jacobian(fun,q),{A,g}=normal(J,r),pg=g.map((x,i)=>(q[i]<=1e-12&&x>0)||(q[i]>=1-1e-12&&x<0)?0:x);
        if(cost<1e-25||Math.max(...pg.map(Math.abs))<1e-12){converged=true;reason='gradient_or_residual';break;}
        let accepted=false;
        for(let trial=0;trial<14;trial++){
          const active=q.map((x,i)=>(x<=1e-10&&g[i]>0)||(x>=1-1e-10&&g[i]<0));
          const B=A.map((row,i)=>row.map((x,j)=>(active[i]||active[j])?(i===j?1:0):x+(i===j?lambda*Math.max(A[i][i],1e-6):0)));
          const delta=linearSolve(B,g.map((x,i)=>active[i]?0:-x));if(!delta){lambda*=10;continue;}
          const q2=q.map((x,i)=>clamp(x+delta[i],0,1)),r2=fun(q2),cost2=dot(r2,r2),move=norm(q2.map((x,i)=>x-q[i]));
          if(Number.isFinite(cost2)&&cost2<cost){q=q2;r=r2;const change=cost-cost2;cost=cost2;lambda=Math.max(1e-14,lambda/3);accepted=true;
            if(move<1e-11||change<1e-25){converged=true;reason='step';}break;}
          lambda*=10;
        }
        if(converged)break;
        if(!accepted){const projected=Math.max(...pg.map(Math.abs));converged=projected<1e-8;reason=converged?'stationary':'no_descent';break;}
      }
      const theta=unpack(q),rfinal=residual(theta),sse=dot(rfinal,rfinal),boundaryHit=q.some(x=>x<1e-5||x>1-1e-5);
      const e=integrateTrajectory(c,theta,{collect:false}).prediction;
      const validEvent=!['invalid_input','integration_failed','outside_model_domain'].includes(e.status)&&!(e.status==='ground_before_goal'&&e.eventTime<times[1]-1e-9);
      results.push({theta,sse,rmse:Math.sqrt(sse/6),converged,accepted:converged&&validEvent,iterations,boundaryHit,reason:validEvent?reason:'ground_before_observation',initialTheta:st.slice(),q});
    }
    const candidates=results.filter(x=>x.accepted).sort((a,b)=>a.sse-b.sse),best=candidates[0];
    if(!best)return {ok:false,theta:null,rmse:null,residuals:[],starts:results.map(({q,...x})=>x),boundaryHit:results.some(x=>x.boundaryHit),diagnostics:{algorithm:'scaled bounded LM',message:'所有初值均未得到有效拟合',freeParameters:indices.map(j=>parameterNames[j])}};
    const r=residual(best.theta),J=jacobian(fun,best.q),{A}=normal(J,r),sv=eigenSym(A).map(x=>Math.sqrt(Math.max(0,x.value))),rank=sv.filter(x=>x>sv[0]*1e-8).length;
    return {ok:true,theta:best.theta,rmse:best.rmse,sse:best.sse,residuals:r,starts:results.map(({q,...x})=>x),boundaryHit:best.boundaryHit,
      diagnostics:{algorithm:'scaled bounded Levenberg–Marquardt, finite differences',termination:best.reason,freeParameters:indices.map(j=>parameterNames[j]),lowerBounds:indices.map(j=>LB[j]),upperBounds:indices.map(j=>UB[j]),singularValues:sv,rank,condition:sv[n-1]>0?sv[0]/sv[n-1]:null,
        parameterScales:indices.map(j=>scales[j]),jacobian:J,rankThreshold:1e-8,weightedResidualSigma:c.sigma>0?c.sigma:null,
        message:'多初值仅核查局部解；零残差不确定测量噪声，b不等于真实转速。'}};
  }
  function fitObservations(input) {
    const c=normalize(input),fit=fitObservationsRaw(c);
    if(fit.ok&&fit.residuals.length){const maxResidual=Math.max(...fit.residuals.map(Math.abs)),threshold=c.sigma>0?5*c.sigma:1e-6;
      fit.diagnostics.compatibility={flag:maxResidual>threshold?(c.sigma>0?'incompatible_with_declared_noise':'large_residual'):'not_flagged_by_residual_rule',maxAbsResidual_m:maxResidual,threshold_m:threshold,
        rule:c.sigma>0?'max absolute coordinate residual > 5 * declared sigma':'max absolute coordinate residual > 1e-6 m; engineering threshold with no noise scale',isFormalStatisticalTest:false};}
    return fit;
  }
  function solve(input) {
    const started=Date.now(),c=normalize(input),issues=validateConfig(c);
    const base={version:VERSION,units:'SI',coordinates:{handedness:'right',goalAxis:'x',upAxis:'z'},config:c,theta:null,trajectory:[],prediction:blankPrediction('invalid_input'),fit:{ok:false,rmse:null,residuals:[],starts:[],boundaryHit:false,diagnostics:{issues}},issues,
      solver:{fit:'classical RK4',fitStep_s:c.fitH,prediction:'classical RK4 with step-local event bisection',predictionStep_s:c.h,eventTimeTolerance_s:1e-12,horizon_s:c.horizon},
      parameterRoles:{free:c.mode==='inverse'?(c.model==='M0'?['vx','vy','vz']:c.model==='M1'?['vx','vy','vz','ax','ay','az']:activeIndices(c).map(j=>['vx','vy','vz','Cd','by','bz'][j])):[],
        fixed:c.mode==='inverse'?['m','R','rho','g','wind','initial','bx=0']:['velocity','Cd','b','constantAcceleration','m','R','rho','g','wind','initial'],
        inactive:[...ADVANCED_FIELDS,...(c.model==='M0'?['Cd','b','wind','rho','constantAcceleration','forceMask']:c.model==='M1'?['Cd','b','wind','rho','g (actual acceleration explicitly prescribed)','forceMask']:c.model==='M2'?['b','constantAcceleration','forceMask']:['constantAcceleration',...(!forceFlags(c).drag?['Cd']:[]),...(!forceFlags(c).magnus?['b']:[])])]}};
    if(issues.length)return base;
    let th=thetaFromConfig(c),fit={ok:true,rmse:null,residuals:[],starts:[],boundaryHit:false,diagnostics:{algorithm:'forward',freeParameters:[],message:'给定参数正向仿真，不使用观测估计。'}};
    if(c.mode==='inverse'){fit=fitObservations(c);if(!fit.ok)return {...base,fit,prediction:blankPrediction('fit_failed'),elapsedMs:Date.now()-started};th=fit.theta;if(fit.constantAcceleration)c.constantAcceleration=fit.constantAcceleration;}
    th=effectiveTheta(c,th);
    if(c.model==='M3'){c.Cd=th[3];c.b=[forceFlags(c).magnus?(c.mode==='inverse'?0:c.b[0]):0,th[4],th[5]];}
    if(fit.diagnostics.compatibility&&['incompatible_with_declared_noise','large_residual'].includes(fit.diagnostics.compatibility.flag))
      issues.push('模型与观测存在显著残差：此处仅为有边界的条件性最小二乘输出，不能作为对这些观测的可接受解释；残差警示不是正式统计检验。');
    const result=integrateTrajectory(c,th);
    const predictedObservations=predictPositions(c,th,c.observations.map(o=>o.t)).map((p,i)=>
      ['ground_before_goal','integration_failed','outside_model_domain'].includes(result.prediction.status)&&c.observations[i].t>(result.prediction.eventTime||0)+1e-10?null:p);
    return {...base,config:c,theta:th,...result,fit,predictedObservations,elapsedMs:Date.now()-started,
      assumptions:['刚性球首飞，首次触地终止','静态均匀风；质量、半径与空气密度为配置假设','M3的b为有效参数；反演固定bx=0','数值验证不等于独立现实实验确认'],
      equation:c.model==='M1'?'r′=v, v′=constantAcceleration':c.model==='M0'?'r′=v, v′=g':c.model==='M2'?'r′=v, v′=g−rho*A*Cd*|w|*w/(2*m)':
        'r′=v, v′=g'+(forceFlags(c).drag?'−rho*A*Cd*|w|*w/(2*m)':'')+(forceFlags(c).magnus?'+b×w':'')};
  }
  function compare(input,comparison='refit') {
    const baseRun=solve(input),runs=[];
    for(const model of ['M0','M1','M2','M3']){const c=normalize(input);c.model=model;
      if(comparison==='fixed'&&baseRun.theta){c.mode='forward';c.velocity=baseRun.theta.slice(0,3);c.Cd=baseRun.theta[3];c.b=[0,...baseRun.theta.slice(4)];}
      else c.mode='inverse';runs.push(solve(c));}
    return {comparison,runs};
  }
  function rng(seed) {let a=seed>>>0;return ()=>{a|=0;a=a+0x6D2B79F5|0;let t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return ((t^t>>>14)>>>0)/4294967296;};}
  function normalRandom(random) {return Math.sqrt(-2*Math.log(Math.max(1e-15,random())))*Math.cos(2*Math.PI*random());}
  function quantiles(xs) {if(!xs.length)return null;const a=xs.slice().sort((a,b)=>a-b);return [.025,.5,.975].map(p=>{const x=p*(a.length-1),i=Math.floor(x);return a[i]+(a[Math.min(i+1,a.length-1)]-a[i])*(x-i);});}
  function describeSamples(samples,method) {
    const statusCounts={crossed_goal_plane:0,ground_before_goal:0,no_crossing_within_horizon:0,invalid_input:0,fit_failed:0,integration_failed:0,outside_model_domain:0,cancelled:0};
    for(const s of samples)statusCounts[s.status]=(statusCounts[s.status]||0)+1;
    const arrivals=samples.filter(s=>s.status==='crossed_goal_plane'),n=arrivals.length,N=samples.length;
    const centerInside=arrivals.filter(s=>s.centerInside).length,clearanceInside=arrivals.filter(s=>s.clearanceInside).length;
    let conditionalEllipse=null;
    if(method==='bootstrap'&&n>=3){const mean=[arrivals.reduce((s,a)=>s+a.yGoal,0)/n,arrivals.reduce((s,a)=>s+a.zGoal,0)/n];
      const covariance=[[0,0],[0,0]];for(const a of arrivals){const d=[a.yGoal-mean[0],a.zGoal-mean[1]];for(let i=0;i<2;i++)for(let j=0;j<2;j++)covariance[i][j]+=d[i]*d[j]/(n-1);}
      const ev=eigenSym(covariance);if(ev[1].value>1e-14)conditionalEllipse={mean,covariance,axes:ev.map(x=>Math.sqrt(5.991464547*x.value)),angle:Math.atan2(ev[0].vector[1],ev[0].vector[0]),chiSquareScale:5.991464547,
        interpretation:'到达样本条件下的bootstrap经验协方差描述椭圆（高斯参考尺度5.991）；不是后验95%区域，未验证覆盖率',denominator:n};}
    return {total:N,completed:N,reached:n,...statusCounts,centerInside,clearanceInside,boundaryHit:samples.filter(s=>s.boundaryHit).length,
      arrivalFraction:N?n/N:null,centerInsideFractionOfArrivals:n?centerInside/n:null,clearanceInsideFractionOfArrivals:n?clearanceInside/n:null,
      denominators:{arrivalFraction:N,centerInsideFractionOfArrivals:n,clearanceInsideFractionOfArrivals:n,conditionalEllipse:n},
      arrivalTimeQuantiles:method==='bootstrap'?quantiles(arrivals.map(s=>s.tGoal)):null,
      yQuantiles:method==='bootstrap'?quantiles(arrivals.map(s=>s.yGoal)):null,zQuantiles:method==='bootstrap'?quantiles(arrivals.map(s=>s.zGoal)):null,
      conditionalEllipse,interpretation:method==='bootstrap'?'围绕拟合模型生成独立高斯坐标噪声并重拟合；频率是该算法与噪声假设下的稳定性，不是现实进球概率':'逐坐标正负扰动的12个确定性情景；频率和范围不具现实概率含义'};
  }
  function computeSample(input,index,seed=20261107,method='bootstrap',baseRun) {
    const c=normalize(input);c.mode='inverse';const base=baseRun||solve(c);
    if(!base.theta)return {index,status:'fit_failed',boundaryHit:false,yGoal:null,zGoal:null,tGoal:null};
    const random=rng((seed+Math.imul(index+1,0x9e3779b1))>>>0);
    c.observations=c.observations.map((o,i)=>({t:o.t,position:(method==='bootstrap'?base.predictedObservations[i]:o.position).map((x,j)=>x+(method==='bootstrap'?c.sigma*normalRandom(random):Math.floor(index/2)===i*3+j?(index%2?1:-1)*c.sigma:0))}));
    c.fitStart=base.theta;const run=solve(c);
    return {index,status:run.prediction.status,...run.prediction,theta:run.theta,observations:c.observations,boundaryHit:run.fit.boundaryHit,
      rmse:run.fit.rmse,fitDiagnostics:run.fit.diagnostics,seed:(seed+Math.imul(index+1,0x9e3779b1))>>>0};
  }
  function uncertainty(input,n=50,seed=20261107,onProgress,method='bootstrap') {
    const c=normalize(input);c.mode='inverse';if(!['bootstrap','scenario'].includes(method))throw new Error('不支持的统计方法');
    if(method==='bootstrap'&&!(c.sigma>0))throw new Error('bootstrap须明确设置sigma>0；零噪声不输出置信标签');
    if(!Number.isInteger(n)||n<1||n>2000)throw new Error('样本数须为1–2000整数');
    const baseRun=solve(c);if(!baseRun.theta)throw new Error('基准拟合失败，无法抽样');
    const total=method==='scenario'?12:n,samples=[];
    for(let i=0;i<total;i++){const sample=computeSample(c,i,seed,method,baseRun);samples.push(sample);if(onProgress)onProgress({completed:i+1,total,sample});}
    return {method,seed,randomAlgorithm:'mulberry32 plus Box–Muller; independent per-index streams',samples,summary:describeSamples(samples,method),baseRun,
      noiseAssumption:{sigma:c.sigma,distribution:method==='bootstrap'?'independent Gaussian per coordinate; initial position exact':'deterministic +/- sigma per coordinate'},version:VERSION};
  }
  return {VERSION,defaultConfig,validateConfig,evaluateForces,integrateTrajectory,predictPositions,fitObservations,solve,run:solve,compare,uncertainty,computeSample,describeSamples};
};
const createDensity=function createDensity() {
  const VERSION = 'upc-density-1.1.0';
  const SQRT_2PI = Math.sqrt(2 * Math.PI);
  const finite = Number.isFinite;
  const failure = (message) => ({version:VERSION,status:'invalid_input',issues:[message]});
  function plainObject(value) { return value !== null && typeof value === 'object' && !Array.isArray(value); }

  // For |x|<3, integrate the Gaussian with a convergent positive-term series.
  // For x>=3, use the Laplace continued fraction for Q(x)/phi(x).
  // No fitted third-party polynomial coefficients are embedded. Returning tails
  // directly avoids 1-CDF cancellation in small, far-tail rectangles.
  function positiveNormalTail(x) {
    if (x === Infinity || x > 38.5) return 0;
    const phi = Math.exp(-0.5*x*x)/SQRT_2PI;
    if (x >= 3) {
      let fraction = 0;
      for (let k = 128; k >= 1; k--) fraction = k/(x+fraction);
      return phi/(x+fraction);
    }
    let term = x, sum = x;
    for (let k = 1; k < 256; k++) {
      term *= x*x/(2*k+1);
      sum += term;
      if (term <= Number.EPSILON * Math.abs(sum)) break;
    }
    return Math.max(0, Math.min(0.5, 0.5-phi*sum));
  }
  function normalCDF(x) {
    if (typeof x !== 'number' || Number.isNaN(x)) return NaN;
    return x >= 0 ? 1-positiveNormalTail(x) : positiveNormalTail(-x);
  }
  function normalInterval(low, high) {
    if (low >= high) return 0;
    let result;
    if (low >= 0) result = positiveNormalTail(low)-positiveNormalTail(high);
    else if (high <= 0) result = positiveNormalTail(-high)-positiveNormalTail(-low);
    else result = 1-positiveNormalTail(high)-positiveNormalTail(-low);
    return Math.max(0, Math.min(1, result));
  }
  function readBounds(bounds, requireFinite) {
    if (!plainObject(bounds)) return null;
    const values = ['yMin','yMax','zMin','zMax'].map(key=>bounds[key]);
    if (values.some(v=>typeof v!=='number'||Number.isNaN(v)||(requireFinite&&!finite(v)))) return null;
    const [yMin,yMax,zMin,zMax] = values;
    if (yMin>yMax||zMin>zMax) return null;
    return {yMin,yMax,zMin,zMax};
  }
  function validateModel(model) {
    if (!plainObject(model)) return 'A density model object is required.';
    if (model.status !== 'ok' && model.status !== 'degenerate') return 'Density model is not usable; inspect its build status.';
    if (model.schema !== 'isotropic-gaussian-kde-1' || !Array.isArray(model.centers) ||
        !Array.isArray(model.weights) || model.weights.length!==model.centers.length ||
        !Array.isArray(model.issues)||model.issues.some(issue=>typeof issue!=='string')) return 'Malformed or unsupported KDE model.';
    if (model.centers.some(p=>!plainObject(p)||!finite(p.y)||!finite(p.z)) || model.weights.some(w=>!finite(w)||w<0)) return 'Invalid centers or weights.';
    let weightSum=0,compensation=0;
    for(const weight of model.weights){const corrected=weight-compensation,next=weightSum+corrected;compensation=(next-weightSum)-corrected;weightSum=next;}
    if (Math.abs(weightSum-(model.centers.length?1:0))>1e-12) return 'Conditional KDE weights must sum to one.';
    if (model.planeArrivalMass !== null) {
      const a = model.planeArrivalMass;
      if (!plainObject(a)||!Number.isSafeInteger(a.nTotal)||a.nTotal<1||a.nArrival!==model.centers.length||a.nArrival>a.nTotal||
          !finite(a.fraction)||Math.abs(a.fraction-a.nArrival/a.nTotal)>1e-15) return 'Invalid all-trial arrival accounting.';
    }
    if(model.status==='degenerate') return model.h===null?null:'Degenerate models must not claim a usable smoothing bandwidth.';
    if(model.centers.length<2||!finite(model.h)||model.h<=0||!finite(2*Math.PI*model.h*model.h)||!finite(1/(2*Math.PI*model.h*model.h))) return 'Malformed or unsupported KDE bandwidth.';
    if(model.supportZMin!==null&&!finite(model.supportZMin)) return 'supportZMin must be null or finite metres.';
    if(!Array.isArray(model.kernelNormalizers)||model.kernelNormalizers.length!==model.centers.length||model.kernelNormalizers.some(n=>!finite(n)||n<=0||n>1)) return 'Invalid per-kernel support normalizers.';
    return null;
  }
  function buildDensity(points, options = {}) {
    if (!Array.isArray(points)||!plainObject(options)) return failure('points must be an array and options an object.');
    if (points.some(p=>!plainObject(p)||!finite(p.y)||!finite(p.z))) return failure('Every point must contain finite y,z in metres; invalid samples are not silently dropped.');
    const allowed = ['sigma_h','minBandwidth','planeArrivalMass','supportZMin'];
    const unknown = Object.keys(options).filter(key=>!allowed.includes(key));
    if (unknown.length) return failure('Unknown options: '+unknown.join(', '));
    const minBandwidth = options.minBandwidth === undefined ? 0.001 : options.minBandwidth;
    const supportZMin = options.supportZMin === undefined ? null : options.supportZMin;
    if(supportZMin!==null&&!finite(supportZMin)) return failure('supportZMin must be null or a finite height in metres.');
    const manual = options.sigma_h !== undefined && options.sigma_h !== null;
    if (!finite(minBandwidth)||minBandwidth<=0 || (manual&&(!finite(options.sigma_h)||options.sigma_h<=0))) return failure('sigma_h and minBandwidth must be positive finite lengths in metres.');
    const centers = points.map(p=>({y:p.y,z:p.z})), n = centers.length;
    let arrival = null;
    if (options.planeArrivalMass !== undefined && options.planeArrivalMass !== null) {
      const a = options.planeArrivalMass;
      if (!plainObject(a)||!Number.isSafeInteger(a.nArrival)||!Number.isSafeInteger(a.nTotal)||
          a.nTotal<1||a.nArrival<0||a.nArrival!==n||a.nArrival>a.nTotal) return failure('planeArrivalMass needs nArrival=points.length and nTotal>=nArrival, with a positive all-trial integer denominator.');
      arrival = {nArrival:a.nArrival,nTotal:a.nTotal,fraction:a.nArrival/a.nTotal,
        interpretation:'Observed plane-arrival fraction among every attempted trial, including fit/integration failures. It is not true arrival or goal probability.'};
    }
    const result = {version:VERSION,schema:'isotropic-gaussian-kde-1',status:'degenerate',issues:[],sampleCount:n,
      centers,weights:n?centers.map(()=>1/n):[],h:null,supportZMin,kernelNormalizers:null,
      units:{coordinates:'m',h:'m',density:'m^-2',integratedMass:'dimensionless'},
      planeArrivalMass:arrival,
      bandwidth:{method:manual?'user_isotropic':'Scott_2D_isotropic',requested:manual?options.sigma_h:null,
        minBandwidth,raw:null,floorApplied:false,formula:'n^(-1/6) * sqrt((sampleVarianceY + sampleVarianceZ)/2)',
        interpretation:'Bandwidth is only a smoothing scale; it is not measurement noise or additional physical uncertainty.'},
      normalization:{support:supportZMin===null?'R^2':'R x [supportZMin, +Infinity)',conditionalIntegral:null,observedMassIntegral:null,truncatedToGoal:false,
        kernelRule:supportZMin===null?'Untruncated kernels':'Each kernel is truncated at supportZMin and separately renormalized; original mixture weights and bandwidth are unchanged.'},
      interpretation:{conditional:'Estimated density conditional on the supplied valid plane-arrival samples.',
        observedMass:'Conditional density times the observed arrival fraction; unavailable without all-trial counts.',
        warning:'Simulation/algorithm-conditional estimate, not calibrated real-world probability or a posterior.'}};
    if (n<2) {
      result.issues.push('At least two arrival samples are required. A zero/single-sample outcome is not an uncertainty density.');
      return result;
    }
    let meanY=0,meanZ=0,m2Y=0,m2Z=0;
    centers.forEach((p,i)=>{const dy=p.y-meanY,dz=p.z-meanZ;meanY+=dy/(i+1);meanZ+=dz/(i+1);m2Y+=dy*(p.y-meanY);m2Z+=dz*(p.z-meanZ);});
    const varianceY=m2Y/(n-1), varianceZ=m2Z/(n-1);
    result.sampleMoments={meanY,meanZ,varianceY,varianceZ,varianceDenominator:n-1};
    if (![meanY,meanZ,varianceY,varianceZ].every(finite) || varianceY<0||varianceZ<0) return failure('Sample moments exceed numerical range.');
    if (centers.every(p=>p.y===centers[0].y&&p.z===centers[0].z)) {
      result.issues.push('All samples coincide. Even a chosen smoothing width must not be displayed as inferred uncertainty.');
      return result;
    }
    const raw = manual ? options.sigma_h : Math.pow(n,-1/6)*Math.sqrt(varianceY/2+varianceZ/2);
    const h = Math.max(raw,minBandwidth);
    if (!finite(h)||h<=0||!finite(2*Math.PI*h*h)||!finite(1/(2*Math.PI*h*h))) return failure('Bandwidth is outside the representable numeric range.');
    const kernelNormalizers=centers.map(p=>supportZMin===null?1:normalInterval((supportZMin-p.z)/h,Infinity));
    if(kernelNormalizers.some(v=>v<=0||!finite(v)||!finite(1/(2*Math.PI*h*h*v)))) return failure('The retained kernel support mass is outside the representable numeric range.');
    result.status='ok';result.h=h;result.bandwidth.raw=raw;result.bandwidth.floorApplied=raw<minBandwidth;
    result.kernelNormalizers=kernelNormalizers;
    result.normalization.conditionalIntegral=1;
    result.normalization.observedMassIntegral=arrival===null?null:arrival.fraction;
    if (raw<minBandwidth) result.issues.push('The declared numerical bandwidth floor was applied; it adds smoothing, not evidence for uncertainty.');
    return result;
  }
  function densityValue(model,y,z) {
    if(model.supportZMin!==null&&z<model.supportZMin)return 0;
    let sum=0;
    for(let i=0;i<model.centers.length;i++) {
      const dy=(y-model.centers[i].y)/model.h,dz=(z-model.centers[i].z)/model.h;
      sum+=model.weights[i]*Math.exp(-0.5*(dy*dy+dz*dz))/model.kernelNormalizers[i];
    }
    return sum/(2*Math.PI*model.h*model.h);
  }
  function evaluateDensity(model,y,z) {
    const error=validateModel(model);if(error)return failure(error);
    if (!finite(y)||!finite(z)) return failure('Evaluation coordinates must be finite metres.');
    if(model.status==='degenerate')return {status:'degenerate',conditionalDensity:null,observedMassDensity:null,issues:model.issues.slice()};
    const value=densityValue(model,y,z);
    return {status:'ok',conditionalDensity:value,observedMassDensity:model.planeArrivalMass===null?null:value*model.planeArrivalMass.fraction,units:'m^-2'};
  }
  function integrateRectangle(model,bounds) {
    const error=validateModel(model);if(error)return failure(error);
    const b=readBounds(bounds,false);if(!b)return failure('Rectangle needs ordered numeric yMin,yMax,zMin,zMax; infinities are allowed.');
    if(model.status==='degenerate')return {status:'degenerate',conditionalMass:null,observedMass:null,issues:model.issues.slice()};
    let value=0;
    const zLower=model.supportZMin===null?b.zMin:Math.max(b.zMin,model.supportZMin);
    for(let i=0;i<model.centers.length;i++) {
      const p=model.centers[i];
      value+=model.weights[i]*normalInterval((b.yMin-p.y)/model.h,(b.yMax-p.y)/model.h)*normalInterval((zLower-p.z)/model.h,(b.zMax-p.z)/model.h)/model.kernelNormalizers[i];
    }
    value=Math.max(0,Math.min(1,value));
    return {status:'ok',bounds:b,conditionalMass:value,observedMass:model.planeArrivalMass===null?null:value*model.planeArrivalMass.fraction,
      units:'dimensionless',method:'Separable normal CDF intervals, respecting per-kernel supportZMin normalization; no renormalization to the queried rectangle.'};
  }
  function gridDensity(model,bounds,nx=61,ny=61) {
    const error=validateModel(model);if(error)return failure(error);
    const b=readBounds(bounds,true);
    if(!b||b.yMin===b.yMax||b.zMin===b.zMax||!finite(b.yMax-b.yMin)||!finite(b.zMax-b.zMin)) return failure('Grid requires finite, nonzero rectangle dimensions.');
    if(!Number.isSafeInteger(nx)||!Number.isSafeInteger(ny)||nx<2||ny<2||nx*ny>1000000) return failure('Grid dimensions must be integers >=2, with at most 1,000,000 points.');
    if(model.status==='degenerate')return {status:'degenerate',conditionalDensity:null,observedMassDensity:null,issues:model.issues.slice()};
    const y=Array.from({length:nx},(_,i)=>b.yMin+(b.yMax-b.yMin)*(i/(nx-1)));
    const z=Array.from({length:ny},(_,i)=>b.zMin+(b.zMax-b.zMin)*(i/(ny-1)));
    const conditionalDensity=z.map(zv=>y.map(yv=>densityValue(model,yv,zv)));
    return {status:'ok',bounds:b,nx,ny,y,z,layout:'density[zIndex][yIndex]',conditionalDensity,
      observedMassDensity:model.planeArrivalMass===null?null:conditionalDensity.map(row=>row.map(v=>v*model.planeArrivalMass.fraction)),
      rectangleMass:integrateRectangle(model,b),units:{coordinates:'m',density:'m^-2'},
      note:'Finite plotting window only; values are not renormalized to the grid or goal opening.'};
  }
  return {version:VERSION,buildDensity,evaluateDensity,integrateRectangle,gridDensity,normalCDF};
};
const createStudy=function createStudy(P,D){
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
};
const createMapping=function createMapping(P,D,S){
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
};
const config={
  "version": "upc-combination-study-0.2",
  "base": {
    "schemaVersion": "1.0",
    "model": "M3",
    "mode": "inverse",
    "initial": [
      0,
      0,
      0.11
    ],
    "observations": [
      {
        "t": 0.05,
        "position": [
          1.5,
          0.45,
          0.44
        ]
      },
      {
        "t": 0.1,
        "position": [
          2.96,
          0.86,
          0.74
        ]
      }
    ],
    "m": 0.43,
    "R": 0.11,
    "rho": 1.225,
    "g": 9.80665,
    "wind": [
      0,
      0,
      0
    ],
    "sigma": 0.01,
    "h": 0.001,
    "fitH": 0.005,
    "horizon": 2,
    "goalX": 11,
    "goalHalfWidth": 3.66,
    "goalHeight": 2.44,
    "velocity": [
      30.4,
      9.4,
      6.9
    ],
    "Cd": 0.3,
    "b": [
      0,
      0,
      0
    ],
    "constantAcceleration": [
      -16,
      -16,
      -12
    ],
    "forceMask": {
      "drag": true,
      "magnus": true
    }
  },
  "sigma": 0.01,
  "sampleCount": 60,
  "assessmentCount": 60,
  "assessmentEnabled": true,
  "seed": 20261107,
  "wind": [
    0,
    2,
    0
  ],
  "schemes": [
    {
      "id": "M0",
      "preset": "gravity",
      "label": "M0 · 重力",
      "forceMask": {
        "drag": false,
        "magnus": false
      },
      "usesWind": false
    },
    {
      "id": "M1",
      "preset": "constantAccel",
      "label": "M1 · 常加速度代数对照",
      "forceMask": {
        "drag": false,
        "magnus": false
      },
      "usesWind": false
    },
    {
      "id": "M2",
      "preset": "drag",
      "label": "M2 · 重力＋阻力",
      "forceMask": {
        "drag": true,
        "magnus": false
      },
      "usesWind": false
    },
    {
      "id": "M3",
      "preset": "magnus",
      "label": "M3 · 重力＋阻力＋Magnus",
      "forceMask": {
        "drag": true,
        "magnus": true
      },
      "usesWind": false
    },
    {
      "id": "G_M",
      "preset": "fields",
      "label": "重力＋Magnus（无阻力）",
      "forceMask": {
        "drag": false,
        "magnus": true
      },
      "usesWind": false
    }
  ],
  "bandwidth": {
    "sigma_h": null,
    "minBandwidth": 0.001
  },
  "selection": {
    "maxNormalizedRMSE": 2,
    "minArrivalFraction": 0.9,
    "maxNumericalFailureFraction": 0.05,
    "minSamples": 20
  },
  "intervalReplicates": 200
};
const P=createPhysics(),D=createDensity(),S=createStudy(P,D),M=createMapping(P,D,S);
function predictFromFrames(p05,p10,patch={}){return M.predictFromFrames(p05,p10,{...config,...patch,base:{...config.base,...patch.base}});}
function selectionSensitivity(p05,p10,options={}){return M.selectionSensitivity(p05,p10,config,options);}
const F=predictFromFrames;
const exportedMapping={F,predictFromFrames,selectionSensitivity,config,version:M.VERSION};
if(typeof module!=='undefined'&&module.exports)module.exports=exportedMapping;
if(typeof globalThis!=='undefined')globalThis.UPCPredictionMapping=exportedMapping;
