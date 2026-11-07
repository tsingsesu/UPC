function report=scenario_checks()
%SCENARIO_CHECKS 对当前同一受力模型检查积分、事件和物理收支。
p=penalty_scenario('defaults');a=penalty_scenario('start',p);
assert(a.fit.feasible);report.observationMaxError=a.fit.maxResidual;
assert(norm(penalty_scenario('sample',a,0).wind-[0;p.windY;0])<1e-12);
q=p;q.windY=-8;b=penalty_scenario('start',q);
report.initialVelocityChange=norm(a.theta(1:3)-b.theta(1:3));assert(report.initialVelocityChange>1e-5);
q.density=1.1;q.Cd=b.params.density*b.params.Cd/q.density;c=penalty_scenario('start',q);
report.densityDragProductDifference=norm(c.endpoint-b.endpoint);assert(report.densityDragProductDifference<1e-6);
steps=[.004 .002 .001];points=zeros(3,3);
for k=1:3
 z=a;z.params.maxStep=steps(k);z=penalty_inverse('propagate',z);points(k,:)=z.endpoint(1:3)';
end
report.steps=steps;report.maxStepDifference=max(vecnorm(points-points(end,:),2,2));assert(report.maxStepDifference<1e-6);
s0=[a.cfg.r0;a.theta(1:3)';0];
sol=ode45(@(t,s)energyRhs(t,s,a),[0 a.tEnd],s0,odeset('RelTol',1e-11,'AbsTol',1e-12,'MaxStep',.001));
E=.5*a.cfg.m*sum(sol.y(4:6,:).^2,1)-a.cfg.m*a.cfg.g(3)*sol.y(3,:);
report.energyWorkResidual=max(abs(E-E(1)-sol.y(7,:)));assert(report.energyWorkResidual<1e-7);
f=a;f.theta(1:3)=[10 0 1];z=penalty_inverse('propagate',f);assert(strcmp(z.status,'ground')&&abs(z.endpoint(3)-.11)<1e-8);
report.groundEventHeight=z.endpoint(3);
f=a;f.theta(1:3)=[80 0 7];z=penalty_inverse('propagate',f);assert(~z.admissibility.accepted);
f=a;f.theta(5:6)=[1.4 1.4];z=penalty_inverse('propagate',f);assert(~z.admissibility.accepted);
report.implausibleSpeedAndForceRejected=true;
% 用相同非线性方程的中心差分检查解析变化方程及穿越时刻修正。
lin=penalty_inverse('linearize',a);Gfd=zeros(2,8);
for j=1:8
 values=[a.theta([1 2 3 5 6]) a.params.windY a.params.density a.params.Cd];
 h=2e-5*max(1,abs(values(j)));plus=a;minus=a;
 if j<=3,plus.theta(j)=plus.theta(j)+h;minus.theta(j)=minus.theta(j)-h;
 elseif j<=5,k=j+1;plus.theta(k)=plus.theta(k)+h;minus.theta(k)=minus.theta(k)-h;
 elseif j==6,plus.cfg.w(2)=plus.cfg.w(2)+h;minus.cfg.w(2)=minus.cfg.w(2)-h;
 elseif j==7
  plus.cfg.rho=plus.cfg.rho+h;minus.cfg.rho=minus.cfg.rho-h;
  plus.cfg.K=plus.cfg.K*plus.cfg.rho/a.cfg.rho;minus.cfg.K=minus.cfg.K*minus.cfg.rho/a.cfg.rho;
 else,plus.theta(4)=plus.theta(4)+h;minus.theta(4)=minus.theta(4)-h;
 end
 plus=penalty_inverse('propagate',plus);minus=penalty_inverse('propagate',minus);
 assert(strcmp(plus.status,'plane')&&strcmp(minus.status,'plane'));
 Gfd(:,j)=(plus.endpoint(2:3)-minus.endpoint(2:3))/(2*h);
end
report.goalSensitivityMaxDifference=max(abs(Gfd-lin.goalJacobian),[],'all');
assert(report.goalSensitivityMaxDifference<1e-5);
report.linearization=lin;
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
path=fullfile(root,'data','派生数据','罚球结果.json');d=jsondecode(fileread(path));d.verification=report;
fid=fopen(path,'w','n','UTF-8');fprintf(fid,'%s',jsonencode(d,'PrettyPrint',true));fclose(fid);result=d;save(fullfile(root,'data','派生数据','罚球结果.mat'),'result');
fid=fopen(fullfile(fileparts(mfilename('fullpath')),'data','实时情景验证.json'),'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));fclose(fid);disp(report);
end
function ds=energyRhs(t,s,a)
[d,f]=ball_rhs(t,s(1:6),a.theta,a.cfg);ds=[d;dot(f.drag+f.magnus,s(4:6))];
end
