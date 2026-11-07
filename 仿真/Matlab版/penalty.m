function result=penalty(mode,obs)
%PENALTY 单一全程受力模型：根据两次位置观测求参数并预测球门位置。
if nargin<1,mode='full';end
if nargin<2,obs=[1.50 .45 .44;2.96 .86 .74];end
assert(any(strcmp(mode,{'full','quick'})),'模式为full或quick。');
validateattributes(obs,{'numeric'},{'size',[2 3],'finite','real'});
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
data=fullfile(root,'data','派生数据');if ~exist(data,'dir'),mkdir(data);end
clock=tic;
f=penalty_family(obs,.35,mode);
p=penalty_inverse('defaults');p.windY=8;p.Cd=.35;
selected=penalty_inverse('propagate',penalty_inverse('fit',p,obs));
assert(selected.fit.feasible,'当前示例未通过观测与物理筛选。');
selected=rmfield(selected,{'sol','cfg'});
result=struct('mode',mode,'software',struct('matlabVersion',version), ...
 'constants',struct('R',.11,'m',.43,'g',[0;0;-9.80665],'rho',1.225,'r0',[0;0;.11]), ...
 'observations',struct('times',[.05 .1],'positions',obs),'selected',selected, ...
 'trajectoryFamily',f,'timing',struct('elapsedSeconds',toc(clock)), ...
 'notes',{{'同一全程受力模型；指定环境后反求初速度和有效侧力。','符合观测的轨迹不带发生概率；范围仅来自列出的搜索设置。','仅进行数值检查，尚无实验确认。'}});
fid=fopen(fullfile(data,'罚球结果.json'),'w','n','UTF-8');assert(fid>=0);
fprintf(fid,'%s',jsonencode(result,'PrettyPrint',true));fclose(fid);
save(fullfile(data,'罚球结果.mat'),'result');
fprintf('[penalty] %d组参数，5mm下%d条符合观测；当前示例(11,%.4f,%.4f)m，t=%.4fs。\n', ...
 f.nCandidates,f.summary(2).nAccepted,selected.endpoint(2),selected.endpoint(3),selected.tEnd);
end
