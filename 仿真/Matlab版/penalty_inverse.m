function out=penalty_inverse(action,varargin)
%PENALTY_INVERSE 全程受力、观测容差约束与环境参数剖面反演。
% 指定(windY,density,Cd)，联合估计(v0x,v0y,v0z,by,bz)。
% 每个剖面保留一个有界最小二乘代表，不声称穷尽可行集。
switch lower(action)
 case 'defaults'
  out=struct('windY',0,'density',1.225,'Cd',.370687007715,'tolerance',.005,'maxStep',.002);
 case 'limits'
  out=struct('speedMin',10,'speedMax',40,'maxTransverseCoefficient',.5, ...
   'massRange',[.410 .450],'circumferenceRange',[.68 .70], ...
   'windRange',[-15 15],'densityRange',[.90 1.40],'CdRange',[.08 .80], ...
   'classification','球规格为官方要求；速度/气动上限为保守模型筛选，非人体或流体绝对上限。', ...
   'sources',{{'https://www.theifab.com/laws/latest/the-ball/', ...
   'https://doi.org/10.1007/BF02844207', ...
   'https://doi.org/10.1080/02640410400021450'}});
 case 'fit'
  p=varargin{1};obs=[1.50 .45 .44;2.96 .86 .74];
  if numel(varargin)>1,obs=varargin{2};end
  validateattributes(obs,{'numeric'},{'size',[2 3],'real','finite'});
  fields={'windY','density','Cd','tolerance','maxStep'};
  for k=1:numel(fields),validateattributes(p.(fields{k}),{'numeric'},{'scalar','real','finite'});end
  assert(p.density>0&&p.Cd>0&&p.tolerance>0&&p.maxStep>0,'非正物理参数/容差。');
  clock=tic;c=configuration(p);lb=[10 -20 0 -1.5 -1.5];ub=[40 20 25 1.5 1.5];
  starts=[30.4 9.4 6.9 -.06 -.36;32 8 7 .1 -.6];
  opts=optimoptions('lsqnonlin','Display','off','FunctionTolerance',1e-11, ...
   'StepTolerance',1e-11,'OptimalityTolerance',1e-9,'MaxIterations',100, ...
   'FiniteDifferenceType','central','TypicalX',[30 10 7 .2 .4]);
  best=inf;audits=repmat(struct('sse',NaN,'exitflag',0),size(starts,1),1);
  for k=1:size(starts,1)
   [q,sse,~,flag]=lsqnonlin(@(x)residual(x,p,c,obs),starts(k,:),lb,ub,opts);
   audits(k)=struct('sse',sse/1e6,'exitflag',flag);
   if isfinite(sse)&&flag>0&&sse<best,best=sse;qBest=q;flagBest=flag;end
  end
  assert(isfinite(best),'penalty_inverse:SolverFailure','全部初值未收敛，不能报告有效候选。');
  theta=[qBest(1:3) p.Cd qBest(4:5)];
  sol=ode45(@(t,s)ball_rhs(t,s,theta,c),[0 .1],[c.r0;theta(1:3)'], ...
   odeset('RelTol',1e-11,'AbsTol',1e-12,'MaxStep',p.maxStep));
  states=deval(sol,[.05 .1]);positions=states(1:3,:)';r=positions-obs;
  maxResidual=max(abs(r(:)));boundary=any(abs(qBest-lb)<1e-5|abs(qBest-ub)<1e-5);
  fit=struct('maxResidual',maxResidual,'rmse',sqrt(mean(r(:).^2)), ...
   'residual',r,'observedPositions',positions,'feasible',maxResidual<=p.tolerance,'observationCompatible',maxResidual<=p.tolerance, ...
   'exitflag',flagBest,'boundaryHit',boundary,'allStarts',audits, ...
   'kineticEnergy',.5*c.m*sum(theta(1:3).^2),'elapsedSeconds',toc(clock), ...
   'note','逐坐标容差筛选；未通过只表示当前剖面代表未通过，不证明全局不可行。');
  out=struct('params',p,'cfg',c,'theta',theta,'fit',fit);
 case 'propagate'
  f=varargin{1};c=f.cfg;th=f.theta;
  opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',f.params.maxStep, ...
   'Events',@(t,s)events(t,s,c.R));
  sol=ode45(@(t,s)ball_rhs(t,s,th,c),[0 2],[c.r0;th(1:3)'],opts);
  assert(all(isfinite(sol.y(:)))&&sol.x(end)>0,'积分失败。');
  status='timeout';if ~isempty(sol.ie),if sol.ie(end)==1,status='plane';else,status='ground';end,end
  ep=deval(sol,sol.x(end));
  out=f;out.sol=sol;out.tEnd=sol.x(end);out.endpoint=ep;out.status=status;
  out.isOpening=strcmp(status,'plane')&&abs(ep(2))<3.66-c.R&&ep(3)>=c.R&&ep(3)<2.44-c.R;
  out.initialSpeed=norm(th(1:3));
  out.admissibility=assess(out);
  out.fit.physicalAdmissible=out.admissibility.accepted;
  out.fit.feasible=out.fit.observationCompatible&&out.admissibility.accepted;
  out.fit.rejectionReasons=out.admissibility.reasons;
  if ~out.fit.observationCompatible,out.fit.rejectionReasons{end+1}='observation_tolerance';end
 otherwise,error('未知反演操作：%s',action);
end
end

function c=configuration(p)
c=struct('R',.11,'m',.43,'rho',p.density,'g',[0;0;-9.80665], ...
 'r0',[0;0;.11],'w',[0;p.windY;0],'bx',0,'cdExponent',0, ...
 'cdRefSpeed',30,'bDecayTau',inf,'K',p.density*pi*.11^2/(2*.43));
end

function r=residual(q,p,c,obs)
% 固定RK4 h=.0025s用于优化；最终接受检查使用ode45。
h=.0025;s=[c.r0;q(1:3)'];pos=zeros(2,3);b=[0;q(4:5)'];d=c.K*p.Cd;
for k=1:40
 k1=fast(s,c,d,b);k2=fast(s+h*k1/2,c,d,b);
 k3=fast(s+h*k2/2,c,d,b);k4=fast(s+h*k3,c,d,b);
 s=s+h*(k1+2*k2+2*k3+k4)/6;
 if k==20,pos(1,:)=s(1:3)';elseif k==40,pos(2,:)=s(1:3)';end
end
r=1000*(pos(:)-obs(:));
end

function ds=fast(s,c,d,b)
u=s(4:6)-c.w;
lift=[b(2)*u(3)-b(3)*u(2);b(3)*u(1);-b(2)*u(1)];
ds=[s(4:6);c.g-d*norm(u)*u+lift];
end

function [v,stop,dir]=events(t,s,R)
v=[s(1)-11;s(3)-R];if t<1e-9,v(2)=max(v(2),1e-8);end
stop=[1;1];dir=[1;-1];
end

function a=assess(s)
L=penalty_inverse('limits');ts=linspace(0,s.tEnd,max(2,ceil(s.tEnd/.001)+1));
vals=deval(s.sol,ts);u=vals(4:6,:)-s.cfg.w;U=sqrt(sum(u.^2,1));
b=repmat([0;s.theta(5:6)'],1,numel(ts));
CL=sqrt(sum(cross(b,u,1).^2,1))./(s.cfg.K*max(U.^2,1e-16));
reasons={};v=norm(s.theta(1:3));
if v<L.speedMin||v>L.speedMax,reasons{end+1}='launch_speed_guard';end
if s.theta(1)<=0||s.theta(3)<0,reasons{end+1}='launch_direction';end
if max(CL)>L.maxTransverseCoefficient,reasons{end+1}='transverse_coefficient_guard';end
if s.cfg.m<L.massRange(1)||s.cfg.m>L.massRange(2)|| ...
  2*pi*s.cfg.R<L.circumferenceRange(1)||2*pi*s.cfg.R>L.circumferenceRange(2)
 reasons{end+1}='ball_specification';
end
p=s.params;
if p.windY<L.windRange(1)||p.windY>L.windRange(2)||p.density<L.densityRange(1)|| ...
 p.density>L.densityRange(2)||p.Cd<L.CdRange(1)||p.Cd>L.CdRange(2)
 reasons{end+1}='outside_declared_environment_domain';
end
if min(vals(3,:))<s.cfg.R-1e-8,reasons{end+1}='below_ground';end
if s.tEnd<.1,reasons{end+1}='terminal_before_last_observation';end
if strcmp(s.status,'timeout'),reasons{end+1}='no_terminal_event';end
% 不把门框范围当作反演约束；射偏也是有效物理结果。
a=struct('accepted',isempty(reasons),'reasons',{reasons}, ...
 'maxTransverseCoefficient',max(CL),'initialSpeed',v, ...
 'initialElevationDegrees',atan2d(s.theta(3),norm(s.theta(1:2))), ...
 'coefficientCheckMaxSpacing',s.tEnd/(numel(ts)-1),'limits',L);
end
