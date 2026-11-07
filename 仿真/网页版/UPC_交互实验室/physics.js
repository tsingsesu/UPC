/* UPC numerical core. Original implementation of the team's M0–M3 equations.
   SI, right handed, x goalward, z up. Rendering never controls integration. */
function createPhysics() {
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
}
if(typeof globalThis!=='undefined') globalThis.UPC=createPhysics();
if(typeof module!=='undefined'&&module.exports) module.exports={createPhysics,...createPhysics()};
