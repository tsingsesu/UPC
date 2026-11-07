function report=scenario_checks()
%SCENARIO_CHECKS 全程反演、可辨识性、能量与事件的数值验证。
p=penalty_scenario('defaults');a=penalty_scenario('start',p);
report.baselineEndpoint=a.endpoint';report.baselineDifference=norm(a.endpoint-a.reference.endpoint);
assert(report.baselineDifference<1e-7&&a.fit.feasible);
p.windY=8;b=penalty_scenario('start',p);
report.windAtLaunch=penalty_scenario('sample',b,0).wind';
assert(report.windAtLaunch(2)==8,'横风必须从t=0参与动力学。');
report.refittedInitialVelocityChange=norm(a.theta(1:3)-b.theta(1:3));
assert(report.refittedInitialVelocityChange>1e-5);
assert(b.fit.feasible&&b.fit.maxResidual<=p.tolerance);
report.windFitResidual=b.fit.residual;
% rho*Cd不变应产生同一有效力场与反演结果，而非双算独立证据。
q=p;q.density=1.1;q.Cd=p.density*p.Cd/q.density;c=penalty_scenario('start',q);
report.densityDragProductDifference=norm(c.endpoint-b.endpoint);
assert(report.densityDragProductDifference<1e-6);
% 固定反演参数，只改变积分设置，避免把重拟合补偿误当积分收敛。
steps=[.004 .002 .001];points=zeros(3,3);
for k=1:3
 f=b;f.params.maxStep=steps(k);s=penalty_inverse('propagate',f);points(k,:)=s.endpoint(1:3)';
end
report.steps=steps;report.stepEndpoints=points;
report.maxStepDifference=max(vecnorm(points-points(end,:),2,2));assert(report.maxStepDifference<1e-6);
[~,forces]=ball_rhs(.03,deval(b.sol,.03),b.theta,b.cfg);u=forces.relativeVelocity;
report.relativeDragPower=dot(forces.drag,u);report.relativeMagnusPower=dot(forces.magnus,u);
assert(report.relativeDragPower<0&&abs(report.relativeMagnusPower)<1e-10);
% 独立累积空气动力对地面系做功，不能强令有风机械能单调下降。
s0=[b.cfg.r0;b.theta(1:3)';0];
sol=ode45(@(t,s)energyRhs(t,s,b),[0 b.tEnd],s0,odeset('RelTol',1e-11,'AbsTol',1e-12,'MaxStep',.001));
E=.5*b.cfg.m*sum(sol.y(4:6,:).^2,1)-b.cfg.m*b.cfg.g(3)*sol.y(3,:);
report.energyWorkResidual=max(abs(E-E(1)-sol.y(7,:)));assert(report.energyWorkResidual<1e-7);
% 错误候选可见但不纳入族；误差统计用六个坐标最大值而非RMSE替代。
p.Cd=.6;bad=penalty_scenario('start',p);report.rejectedResidual=bad.fit.maxResidual;
assert(~bad.fit.feasible);
% 同一批候选的容差筛选必须嵌套。
A=penalty_scenario('family',.002);B=penalty_scenario('family',.005);C=penalty_scenario('family',.01);
assert(numel(A)<=numel(B)&&numel(B)<=numel(C));report.familyCounts=[numel(A) numel(B) numel(C)];
assert(~isempty(B)&&all([B.maxResidual]<=.005));
% 触地与平面定位都必须定位，不能继续画地下轨迹。
f=a;f.theta(1:3)=[10 0 1];s=penalty_inverse('propagate',f);
assert(strcmp(s.status,'ground')&&abs(s.endpoint(3)-.11)<1e-8);report.groundEvent=s.endpoint';
% 新界面不允许把不同环境的两个解拼接成同次射门。
rejected=false;try,penalty_scenario('branch',a,.2,p);catch e,rejected=strcmp(e.identifier,'penalty_scenario:RefitRequired');end
assert(rejected);report.oldBranchRejected=true;
f=a;f.theta(1:3)=[80 0 7];unphysical=penalty_inverse('propagate',f);
assert(~unphysical.admissibility.accepted&&any(strcmp(unphysical.admissibility.reasons,'launch_speed_guard')));
f=a;f.theta(5:6)=[1.4 1.4];unphysical=penalty_inverse('propagate',f);
assert(~unphysical.admissibility.accepted&&any(strcmp(unphysical.admissibility.reasons,'transverse_coefficient_guard')));
report.implausibleSpeedAndForceRejected=true;
ranges=[-15 .90 .08;15 1.40 .80;-12 1.0 .2;12 1.3 .5];
report.expandedRangeChecks=struct([]);
for k=1:size(ranges,1)
 q=penalty_scenario('defaults');q.windY=ranges(k,1);q.density=ranges(k,2);q.Cd=ranges(k,3);
 z=penalty_scenario('start',q);
 assert(~any(strcmp(z.admissibility.reasons,'outside_declared_environment_domain')));
 assert(z.fit.feasible==(z.fit.observationCompatible&&z.fit.physicalAdmissible));
 report.expandedRangeChecks=[report.expandedRangeChecks;struct('windY',q.windY,'density',q.density,'Cd',q.Cd, ...
 'accepted',z.fit.feasible,'residual',z.fit.maxResidual,'reasons',{z.fit.rejectionReasons})];
end
root=fileparts(mfilename('fullpath'));fid=fopen(fullfile(root,'data','实时情景验证.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));fclose(fid);disp(report);
end
function ds=energyRhs(t,s,b)
[d,f]=ball_rhs(t,s(1:6),b.theta,b.cfg);ds=[d;dot(f.drag+f.magnus,s(4:6))];
end
