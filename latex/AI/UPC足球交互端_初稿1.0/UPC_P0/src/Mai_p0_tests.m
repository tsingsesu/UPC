function report = Mai_p0_tests()
%MAI_P0_TESTS 现有方程的解析、事件、退化和反演回归；不写文件。
started=tic;
cfg=Mai_penalty('api','config');
names={}; errors=[];
theta=[25,4,6,0,0,0]; times=[.05,.1,.2];
sol=Mai_penalty('api','positions',struct('theta',theta,'cfg',cfg,'times',times));
exact=cfg.r0'+times(:)*theta(1:3)+.5*times(:).^2*cfg.g';
e=max(abs(sol.positions-exact),[],'all'); assert(e<1e-10);
names{end+1}='V01 重力解析位置'; errors(end+1)=e;
p=Mai_penalty('api','predict',struct('theta',theta,'cfg',cfg));
exactGoal=cfg.r0'+theta(1:3)*(.44)+.5*cfg.g'*.44^2;
e=max(abs([p.tGoal-.44,p.yGoal-exactGoal(2),p.zGoal-exactGoal(3)]));
assert(e<1e-8&&strcmp(p.status,'goal'));
names{end+1}='V01 起飞与重力过门事件'; errors(end+1)=e;
dragCfg=cfg; dragCfg.g=zeros(3,1); thetaD=[18,0,0,.25,0,0];
s=Mai_penalty('api','positions',struct('theta',thetaD,'cfg',dragCfg,'times',times));
exactX=log(1+cfg.K*.25*18*times)/(cfg.K*.25);
e=max(abs(s.positions(:,1)'-exactX)); assert(e<1e-8);
names{end+1}='V02 一维平方阻力解析'; errors(end+1)=e;
thetaM=[25,4,6,.25,.2,-.3]; state=[cfg.r0;thetaM(1:3)'];
f=Mai_penalty('api','forces',struct('theta',thetaM,'cfg',cfg,'state',state));
e=abs(f.magnus'*state(4:6)); assert(e<1e-10&&f.drag'*state(4:6)<0);
names{end+1}='V03/V04 功率符号与Magnus正交'; errors(end+1)=e;
zeroCfg=cfg; zeroCfg.w=state(4:6); zeroCfg.cdExponent=-.5;
f=Mai_penalty('api','forces',struct('theta',thetaM,'cfg',zeroCfg,'state',state));
e=norm(f.acceleration-cfg.g); assert(e<1e-12&&all(isfinite(f.acceleration)));
names{end+1}='V05 零相对速度含旧幂律情景'; errors(end+1)=e;
early=Mai_penalty('api','predict',struct('theta',[25,0,.1,0,0,0],'cfg',cfg));
e=abs(early.eventTime-2*.1/9.80665);
assert(strcmp(early.status,'ground')&&early.eventTime<.1&&isnan(early.tGoal)&&e<1e-8);
names{end+1}='V06 观测前触地'; errors(end+1)=e;
fallCfg=cfg; fallCfg.r0=[0;0;cfg.R+.01];
fall=Mai_penalty('api','predict',struct('theta',[30,0,-1,0,0,0],'cfg',fallCfg));
fallExact=(-1+sqrt(1+2*9.80665*.01))/9.80665;
e=abs(fall.eventTime-fallExact); assert(strcmp(fall.status,'ground')&&e<1e-9);
names{end+1}='V06 与独立Python共享下降触地解析算例'; errors(end+1)=e;
rest=Mai_penalty('api','predict',struct('theta',[25,0,0,0,0,0],'cfg',cfg));
assert(strcmp(rest.status,'ground')&&rest.eventTime==0&&isnan(rest.tCross));
names{end+1}='V06 地面非起飞初态'; errors(end+1)=0;
timeoutCfg=cfg; timeoutCfg.g=zeros(3,1); timeoutCfg.r0=[0;0;1];
p=Mai_penalty('api','predict',struct('theta',[1,0,0,0,0,0],'cfg',timeoutCfg));
assert(strcmp(p.status,'timeout')&&p.eventTime==2&&isnan(p.tGoal));
names{end+1}='V06 未达门不冒充过门'; errors(end+1)=0;
bad=cfg; bad.R=.12; rejected=false;
try, Mai_penalty('api','config',struct('cfg',bad)); catch, rejected=true; end
assert(rejected); names{end+1}='V06 半径不静默改题给球心';errors(end+1)=0;
p=Mai_penalty('api','analytic',struct('v0',[30.4,9.4,6.9],'a',[-16,-16,-12],'cfg',cfg));
e=max(abs([p.tGoal-.4050083612,p.yGoal-2.4948244142,p.zGoal-1.9203670565]));
assert(e<1e-9); names{end+1}='M1 常加速度代数回归';errors(end+1)=e;
covOk=Mai_penalty('api','covariance',struct('J',eye(6),'cfg',cfg));
covBad=Mai_penalty('api','covariance',struct('J',ones(6),'cfg',cfg));
covEdge=Mai_penalty('api','covariance',struct('J',eye(6),'cfg',cfg,'boundaryHit',true));
assert(covOk.diagnostics.valid&&~covBad.diagnostics.valid&&~covEdge.diagnostics.valid);
assert(all(isnan(covBad.covTheta),'all')&&all(isnan(covEdge.covTheta),'all'));
e=max(abs(covOk.covTheta-cfg.sigma^2*eye(6)),[],'all');assert(e<1e-12);
names{end+1}='局部协方差尺度/秩/边界门控';errors(end+1)=e;
noNoise=cfg;noNoise.sigma=0;
c=Mai_penalty('api','covariance',struct('J',eye(6),'cfg',noNoise));assert(~c.diagnostics.valid);
names{end+1}='未设噪声不生成95%区间';errors(end+1)=0;
obs=[1.5,.45,.44;2.96,.86,.74];
fit=Mai_penalty('api','fit',struct('obs',obs,'cfg',cfg));
e=max(abs([fit.prediction.yGoal-2.5994405492258013,fit.prediction.zGoal-1.9670046154749097,fit.prediction.tGoal-.40050348611563008]));
assert(fit.fit.exitflag>0&&fit.fit.rmse<1e-8&&e<1e-7);
assert(fit.uncertainty.diagnostics.valid);
names{end+1}='官方数据六参数反演回归';errors(end+1)=e;
obs(1,2)=obs(1,2)+.01;
custom=Mai_penalty('api','fit',struct('obs',obs,'cfg',cfg));
assert(custom.fit.exitflag>0&&abs(custom.prediction.yGoal-fit.prediction.yGoal)>1e-3);
names{end+1}='自定义观测真实重拟合';errors(end+1)=custom.fit.rmse;
steps=[.02,.01,.005]; events=zeros(3,3);
for k=1:3
    c=cfg;c.odeMaxStep=steps(k);
    p=Mai_penalty('api','predict',struct('theta',fit.fit.theta,'cfg',c));
    events(k,:)=[p.yGoal,p.zGoal,p.tGoal];
end
e=max(abs(events-events(end,:)),[],'all'); assert(e<1e-5);
names{end+1}='V07 三档事件步长收敛'; errors(end+1)=e;
failedCfg=cfg; failedCfg.freeMask=false(1,6);
failed=Mai_penalty('api','fit',struct('obs',[1.5,.45,.44;2.96,.86,.74],'cfg',failedCfg));
assert(failed.fit.exitflag<=0&&strcmp(failed.prediction.status,'fit_fail'));
assert(~failed.uncertainty.diagnostics.valid&&all(isnan(failed.uncertainty.covTheta),'all'));
names{end+1}='拟合无有效候选返回明确失败'; errors(end+1)=0;
underground=[25,0,-.1,.25,-1,0];
illegalCfg=cfg; illegalCfg.freeMask=[true,false,false,false,false,false];
illegalCfg.thetaFixed=underground;
raw=Mai_penalty('api','positions',struct('theta',underground,'cfg',cfg,'times',[.05,.1]));
assert(all(raw.positions(:,3)>cfg.R));
rejectedFit=Mai_penalty('api','fit',struct('obs',raw.positions,'cfg',illegalCfg));
assert(strcmp(rejectedFit.prediction.status,'fit_fail')&&all([rejectedFit.allStarts.exitflag]==-4));
assert(all([rejectedFit.allStarts.optimizerExitflag]>0));
names{end+1}='数学延拓可拟合但首飞前触地的候选被拒绝';errors(end+1)=0;
noNoiseRun=Mai_penalty('quick',[],struct('writeOutputs',false,'makeFigures',false, ...
    'configOverrides',struct('sigma',0)));
assert(~noNoiseRun.uncertainty.diagnostics.valid&&all(isnan(noNoiseRun.uncertainty.marginal95HalfWidth_yz)));
names{end+1}='完整quick无噪声区间保持缺失而非零';errors(end+1)=0;
boundaryRun=Mai_penalty('quick',[],struct('writeOutputs',false,'makeFigures',false, ...
    'configOverrides',struct('ub',[50,20,25,.2,1.5,1.5])));
assert(boundaryRun.fit.boundaryHit&&~boundaryRun.uncertainty.diagnostics.valid);
assert(all(isnan(boundaryRun.uncertainty.marginal95HalfWidth_yz)));
assert(all(isnan(vertcat(boundaryRun.noiseScenarios.marginal95HalfWidth_yz))));
names{end+1}='完整quick边界区间保持缺失而非零';errors(end+1)=0;
report=struct('passed',true,'testCount',numel(names),'names',{names}, ...
    'measuredErrors',errors,'unitsNote','各测试误差按名称对应位置m/时间s/力功率等；非现实精度', ...
    'elapsedSeconds',toc(started),'matlabVersion',version,'physicalValidation','未完成实验确认');
disp(report);
end
