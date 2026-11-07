function report=scenario_checks()
%SCENARIO_CHECKS 有物理含义的扩展模型验证；不改正式论文结果。
p=penalty_scenario('defaults');a=penalty_scenario('start',p);
report.baselineEndpoint=a.endpoint';report.baselineTime=a.tEnd;
report.baselineDifference=norm(a.endpoint-a.reference.endpoint);
assert(report.baselineDifference<1e-7);
g=p;g.gustY=12;g.gustStart=.15;g.gustDuration=.2;
b=penalty_scenario('start',g);
errors=zeros(1,3);
for k=1:3
    t=[.05,.10,.14];x=penalty_scenario('sample',a,t(k));y=penalty_scenario('sample',b,t(k));
    errors(k)=norm(x.state-y.state);
end
assert(max(errors)<1e-8,'阵风发生前不应改变观测或轨迹。');
report.sameObservationErrors=errors;
report.gustEndpoint=b.endpoint';report.gustShift_m=norm(b.endpoint(2:3)-a.endpoint(2:3));
assert(report.gustShift_m>.01);

q=p;q.windY=-12;c=penalty_scenario('branch',a,.22,q);
report.branchPastError=0;
for t=[0,.05,.1,.15,.21999,.22]
    before=penalty_scenario('sample',a,t);after=penalty_scenario('sample',c,t);
    report.branchPastError=max(report.branchPastError,norm(before.state-after.state));
end
assert(report.branchPastError<1e-10,'改参不应改写过去或跳变状态。');
q.windY=10;cc=penalty_scenario('branch',c,.28,q);
history=penalty_scenario('sample',cc,.24);future=penalty_scenario('sample',cc,.30);
assert(abs(history.wind(2)+12)<1e-12&&abs(future.wind(2)-10)<1e-12);
report.historyWindAndFuture=[history.wind';future.wind'];

state=penalty_scenario('sample',b,.24);u=state.state(4:6)-state.wind;
report.dragRelativePower=dot(state.forces.drag,u);
report.magnusRelativePower=dot(state.forces.magnus,u);
assert(report.dragRelativePower<=0&&abs(report.magnusRelativePower)<1e-10);

last=b.segments(end);cfg=last.cfg;theta=last.theta;s0=deval(last.sol,last.t0);
ext=ode45(@(t,s)energyRhs(t,s,theta,cfg),[last.t0,last.tf],[s0;0], ...
    odeset('RelTol',1e-10,'AbsTol',1e-12,'MaxStep',.001));
energy=.5*cfg.m*sum(ext.y(4:6,:).^2,1)-cfg.m*cfg.g(3)*ext.y(3,:);
report.windWorkBalance_J=max(abs(energy-energy(1)-ext.y(7,:)));
assert(report.windWorkBalance_J<1e-7,'含风时必须核对地面系气动力功。');

steps=[.004,.002,.001];endpoints=zeros(3,4);
for k=1:3
    r=g;r.maxStep=steps(k);sol=penalty_scenario('start',r);
    endpoints(k,:)=[sol.tEnd,sol.endpoint(1:3)'];
end
report.maxSteps=steps;report.refinement=endpoints;
report.refinementDifference_m=max(vecnorm(endpoints(:,2:4)-endpoints(end,2:4),2,2));
assert(report.refinementDifference_m<1e-6);

free=p;free.mode='free';free.speedScale=.6;free.elevationDeg=-10;
ground=penalty_scenario('start',free);
assert(strcmp(ground.status,'ground'),'提前触地须分类且终止。');
report.groundTest=struct('time',ground.tEnd,'point',ground.endpoint(1:3)');
tail=p;tail.windX=10;head=p;head.windX=-10;
ta=penalty_scenario('start',tail);he=penalty_scenario('start',head);
assert(ta.tEnd<a.tEnd&&he.tEnd>a.tEnd);
report.tailCalmHeadTimes=[ta.tEnd,a.tEnd,he.tEnd];

% 一次早期调参在conditioned模式下最早从0.10秒生效。
early=penalty_scenario('branch',a,.03,q);
earlyObs=penalty_scenario('sample',early,.1);original=penalty_scenario('sample',a,.1);
assert(norm(earlyObs.state-original.state)<1e-10);
bad=p;bad.density=nan;rejected=false;
try,penalty_scenario('start',bad);catch,rejected=true;end
assert(rejected,'应拒绝NaN参数。');
root=fileparts(mfilename('fullpath'));path=fullfile(root,'data','实时情景验证.json');
report.passed=true;report.scope='数值和因果性检验；非实测气动标定，非实验确认。';
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));fclose(fid);
disp(report);
end

function ds=energyRhs(t,s,theta,cfg)
[rhs,F]=ball_rhs(t,s(1:6),theta,cfg);
ds=[rhs;dot(F.drag+F.magnus+F.buoyancy,s(4:6))];
end
