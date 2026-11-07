function reference = Mai_web_reference()
%MAI_WEB_REFERENCE 由已冻结MATLAB纯接口生成浏览器M0—M3对照，不改正式结果。
% 从本文件所在包定位输出；只写 data/验证基准/浏览器对照.json。
cfg=Mai_penalty('api','config');
obs=[1.5,.45,.44;2.96,.86,.74];
times=[.05,.1];
reference=struct('schemaVersion','1.0','matlabVersion',version, ...
    'sourceCommit','849f45a8e830d45348239c68cccd9f92dce3986d', ...
    'modelCodeVersion','existing-model-p0-repair-1', ...
    'units','SI','coordinates','right-handed; x goalward, z upward', ...
    'statusMapping',struct('goal','crossed_goal_plane','cross_out','crossed_goal_plane', ...
    'ground','ground_before_goal','timeout','no_crossing_within_horizon', ...
    'fit_fail','fit_failed','invalid_input','invalid_input','integration_failed','integration_failed'), ...
    'scope','数值实现对照；不是独立实验确认。M1为常加速度代数延拓。', ...
    'tolerances',struct('position_m',1e-4,'eventTime_s',1e-5,'force_N',1e-8, ...
    'fitObservation_m',1e-5,'note','优化器参数可有微小差异；优先同参数力、位置和事件比较。'));
cases=cell(1,8);
c0=cfg;c0.freeMask=[true,true,true,false,false,false];c0.thetaFixed=zeros(1,6);
f0=Mai_penalty('api','fit',struct('obs',obs,'cfg',c0));
assert(f0.fit.exitflag>0);
cases{1}=makeCase('inverse_M0','M0','inverse',c0,f0.fit.theta,[],obs,f0.fit);
v1=(4*obs(1,:)-obs(2,:)-3*cfg.r0')/.1;
a1=(obs(2,:)-2*obs(1,:)+cfg.r0')/.05^2;
cases{2}=makeCase('inverse_M1','M1','inverse',cfg,[v1,0,0,0],a1,obs,struct());
c2=cfg;c2.freeMask=[true,true,true,true,false,false];c2.thetaFixed=zeros(1,6);
f2=Mai_penalty('api','fit',struct('obs',obs,'cfg',c2));
assert(f2.fit.exitflag>0);
cases{3}=makeCase('inverse_M2','M2','inverse',c2,f2.fit.theta,[],obs,f2.fit);
f3=Mai_penalty('api','fit',struct('obs',obs,'cfg',cfg));
assert(f3.fit.exitflag>0);
theta=f3.fit.theta;
cases{4}=makeCase('inverse_M3','M3','inverse',cfg,theta,[],obs,f3.fit);
t0=theta;t0(4:6)=0;
cases{5}=makeCase('fixed_initial_M0','M0','forward',cfg,t0,[],obs,struct());
t2=theta;t2(5:6)=0;
cases{6}=makeCase('fixed_initial_M2','M2','forward',cfg,t2,[],obs,struct());
cases{7}=makeCase('fixed_initial_M3','M3','forward',cfg,theta,[],obs,struct());
windCfg=cfg;windCfg.w=[0;3;0];
cases{8}=makeCase('wind_y_plus3_M3','M3','forward',windCfg,theta,[],obs,struct());
reference.cases=[cases{:}];
forceCases=cell(1,4);
state=[cfg.r0;theta(1:3)'];
forceCases{1}=makeForce('nominal_initial',cfg,theta,0,state);
forceCases{2}=makeForce('wind_y_plus3_initial',windCfg,theta,0,state);
zeroCfg=cfg;zeroCfg.w=state(4:6);
forceCases{3}=makeForce('zero_relative_speed',zeroCfg,theta,0,state);
noSpin=theta;noSpin(5:6)=0;
forceCases{4}=makeForce('zero_effective_magnus',cfg,noSpin,0,state);
reference.forceCases=[forceCases{:}];
reference.observationTimes=times;
reference.originalObservations=obs;
targetDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','验证基准');
if ~isfolder(targetDir),mkdir(targetDir);end
target=fullfile(targetDir,'浏览器对照.json');
fid=fopen(target,'w','n','UTF-8');
if fid<0,error('Mai_web_reference:WriteFailed','无法创建浏览器对照数据。');end
written=fprintf(fid,'%s',jsonencode(reference,'PrettyPrint',true));
closed=fclose(fid);
if written<1||closed~=0,error('Mai_web_reference:WriteFailed','浏览器对照数据写入失败。');end
fprintf('MATLAB浏览器对照已生成：8个轨迹算例、4个受力算例。\n');
end

function out=makeCase(id,model,mode,cfg,theta,accel,obs,fit)
sampleTimes=[0,.05,.1,.2,.3];
if strcmp(model,'M1')
    p=Mai_penalty('api','analytic',struct('cfg',cfg,'v0',theta(1:3),'a',accel));
    samplePositions=cfg.r0'+sampleTimes(:)*theta(1:3)+.5*sampleTimes(:).^2*accel;
    positions=samplePositions(2:3,:);
else
    p=Mai_penalty('api','predict',struct('cfg',cfg,'theta',theta,'targetX',11));
    sample=Mai_penalty('api','positions',struct('cfg',cfg,'theta',theta,'times',sampleTimes));
    samplePositions=sample.positions;
    positions=samplePositions(2:3,:);
end
assert(any(strcmp(p.status,{'goal','cross_out'})));
bc=browserConfig(cfg,theta,model,mode,obs);
if ~isempty(accel),bc.constantAcceleration=accel;end
mc=cfg;mc.bDecayTau='infinite';
out=struct('id',id,'model',model,'mode',mode,'config',mc,'browserConfig',bc, ...
    'theta',theta,'constantAcceleration',accel,'prediction',p, ...
    'observationsPredicted',positions,'observationResidual',positions-obs, ...
    'sampleTimes',sampleTimes,'samplePositions',samplePositions,'fitDiagnostics',fit);
end

function bc=browserConfig(cfg,theta,model,mode,obs)
ob(1)=struct('t',.05,'position',obs(1,:));
ob(2)=struct('t',.1,'position',obs(2,:));
bc=struct('model',model,'mode',mode,'initial',cfg.r0','observations',ob, ...
    'm',cfg.m,'R',cfg.R,'rho',cfg.rho,'g',-cfg.g(3),'wind',cfg.w', ...
    'sigma',cfg.sigma,'h',.001,'fitH',cfg.h,'horizon',cfg.maxFlightTime, ...
    'goalX',11,'goalHalfWidth',3.66,'goalHeight',2.44, ...
    'velocity',theta(1:3),'Cd',theta(4),'b',[cfg.bx,theta(5:6)], ...
    'constantAcceleration',[-16,-16,-12]);
end

function out=makeForce(id,cfg,theta,t,state)
force=Mai_penalty('api','forces',struct('cfg',cfg,'theta',theta,'t',t,'state',state));
out=struct('id',id,'browserConfig',browserConfig(cfg,theta,'M3','forward', ...
    [1.5,.45,.44;2.96,.86,.74]),'theta',theta,'t',t,'state',state', ...
    'expected',force);
end
