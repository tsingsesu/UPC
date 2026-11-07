function report = p0_validate(mode)
% 一次调用完成自测、主计算、配置读回及隔离性检查；不写论文。
if nargin<1, mode='full'; end
root=fileparts(fileparts(mfilename('fullpath')));
tests=p0_tests();
r=penalty(mode,[],struct('outputRoot',root));
jsonPath=fullfile(root,r.output.dataDirectory,'Mai罚球结果.json');
loaded=jsondecode(fileread(jsonPath));
cfg=penalty('api','config',struct('cfg',loaded.config));
replay=penalty('api','predict',struct('theta',loaded.fit.theta,'cfg',cfg));
assert(max(abs([replay.yGoal-r.prediction.yGoal,replay.zGoal-r.prediction.zGoal,replay.tGoal-r.prediction.tGoal]))<1e-10);
observations=[1.5,.45,.44;2.96,.86,.74];
custom=zeros(6,3);
for j=1:6
    changed=observations; changed(j)=changed(j)+.001;
    v=penalty('api','fit',struct('obs',changed,'cfg',cfg));
    assert(v.fit.exitflag>0&&v.fit.rmse<1e-7);
    custom(j,:)=[v.prediction.yGoal,v.prediction.zGoal,v.prediction.tGoal];
    assert(norm(custom(j,:)-[r.prediction.yGoal,r.prediction.zGoal,r.prediction.tGoal])>1e-6);
end
before=dir(fullfile(root,'data','验证运行'));
stateBefore=rng;
q=penalty('quick',[],struct('writeOutputs',false,'makeFigures',false));
stateAfter=rng;
after=dir(fullfile(root,'data','验证运行'));
assert(isequal({before.name},{after.name})&&isequal(stateBefore,stateAfter));
assert(~q.output.written&&isempty(q.figures));
assert(r.bootstrap.summary.nSamples==numel(r.bootstrap.samples));
assert(r.synthetic.summary.nSamples==numel(r.synthetic.samples));
assert(all(arrayfun(@(s) ~strcmp(s.status,'ground')||isnan(s.tGoal),r.bootstrap.samples)));
assert(r.verification.rk4GravityMaxError<1e-10&&r.verification.dragAnalyticMaxError<1e-8);
assert(r.verification.magnusPowerPass&&r.verification.dragPowerPass);
report=struct('passed',true,'mode',mode,'tests',tests,'resultPath',strrep(fullfile(r.output.dataDirectory,'Mai罚球结果.json'),'\','/'), ...
    'configurationRoundtrip',true,'sixCoordinatesRecomputed',true,'customPredictions_y_z_t',custom, ...
    'noWriteHasNoFileOrRngSideEffects',true,'matlabVersion',version, ...
    'fullP0ProductComplete',false,'experimentalValidationComplete',false);
target=fullfile(root,r.output.dataDirectory,'P0验收.json');
fid=fopen(target,'w','n','UTF-8'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fprintf('P0_VALIDATION_RESULT=%s\n',target);
end
