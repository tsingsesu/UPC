function report=scenario_ui_checks()
%SCENARIO_UI_CHECKS 用真实回调检验全程约束反演界面（本文件只编写，不自动运行）。
%   检查项：唯一三维轴、四个参数控件、t>0 改参从 0 全程重建、非法输入恢复、
%   不合格候选的红色标示、播放保持、重播、关闭清理 timer。
%   注意：参数编辑会触发全程联合反演，实际运行时间取决于引擎；脚本不等待播放结束。

penalty_view;
f=gcf;
cleanup=onCleanup(@()closeIfLive(f)); %#ok<NASGU>
set(f,'Visible','off');
drawnow;

report=struct(); report.passed=false;

% ---------- 唯一三维轴与控件集合 ----------
ax=findall(f,'Type','axes');
assert(numel(ax)==1,'界面必须只有一个三维轴。');
assert(strcmp(get(ax(1),'Tag'),'Trajectory3D'),'唯一轴 Tag 必须是 Trajectory3D。');
assert(~ismember(get(ax(1),'Tag'),{'Conditioned','Free'}),'不得保留旧的二维/模式轴。');

names={'windY','density','Cd','tolerance'};
for k=1:numel(names)
    nm=names{k};
    sl=control(['Param_' nm]); ed=control(['Edit_' nm]);
    assert(strcmp(get(sl,'Style'),'slider'),['Param_' nm ' 必须是滑块。']);
    assert(strcmp(get(ed,'Style'),'edit'),['Edit_' nm ' 必须是编辑框。']);
end
for tag={'PlayButton','ReplayButton','TimeSlider','ViewSelect','FamilyToggle','ForceToggle','ResetButton'}
    control(tag{1});
end
% 只有四个数值参数，旧参数控件不应存在
for tag={'Param_speedScale','Param_mass','Param_radius','Param_gustY','ModeSelect','WindToggle'}
    assert(isempty(findobj(f,'Tag',tag{1})),['不应再出现旧控件 ' tag{1}]);
end

st=control('StatusText');
famNote=control('FamilyNote');
assert(~isempty(strfind(get(famNote,'String'),'置信区间')), ...
    '族说明必须明确“不是置信区间”。');
report.singleAxes=true;
report.fourParameters=true;

% ---------- t>0 时编辑参数：从 0 全程重建 ----------
sl=control('TimeSlider');
set(sl,'Value',0.5*get(sl,'Max'));
invoke(sl);
assert(getappdata(f,'CurrentTime')>0,'滑块定位后当前时间应大于 0。');
old=getappdata(f,'Scenario');
oldV0=old.theta(1:3);

ed=control('Edit_windY');
newWind=3.5;
set(ed,'String',sprintf('%.2f',newWind));
invoke(ed);
sim=getappdata(f,'Scenario');
p=getappdata(f,'Parameters');
assert(getappdata(f,'CurrentTime')==0,'编辑参数后必须从 t=0 重建并暂停。');
assert(abs(p.windY-newWind)<1e-12,'新横风未写入参数。');
assert(norm(sim.theta(1:3)-oldV0)>1e-9,'全程反演后起脚速度应随之变化（重新拟合）。');
curve=penalty_scenario('curve',sim);
assert(abs(curve.t(1))<1e-9,'重建后轨迹必须从 t=0 开始，不得拼接旧历史。');
assert(abs(curve.t(end)-sim.tEnd)<1e-6,'重建后轨迹终点必须等于 tEnd。');
% 早期观测仍应满足容差（合格候选）；不合格结果本身是合法输出，由红色分支覆盖。
early=curve.t<=min(0.12,sim.tEnd);
assert(any(early),'轨迹长度不足以检查早期段。');
if sim.fit.feasible
    assert(sim.fit.maxResidual<=sim.params.tolerance+1e-12, ...
        '合格的全程反演应满足观测容差。');
    if ~isempty(sim.fit.residual)
        r0=max(abs(sim.fit.residual(1,:)));
        assert(r0<=sim.params.tolerance+1e-12,'最早观测点应满足容差。');
    end
end
assert(sim.fit.feasible,'温和横风候选应满足默认容差。');
assert(norm([get(findobj(ax(1),'Tag','Observed'),'XData')]-[1.5 2.96])<1e-12,'黑色观测点不能随拟合移动。');
report.feasibleAfterEdit=sim.fit.feasible;
report.rebuildFromZero=true;
report.v0Changed=true;
report.earlyWithinTolerance=true;
report.oldHistoryNotAppended=true;

% ---------- 非法输入恢复 ----------
before=getappdata(f,'Parameters');
ed=control('Edit_Cd');
set(ed,'String','NaN');
invoke(ed);
after=getappdata(f,'Parameters');
assert(abs(after.Cd-before.Cd)<1e-12,'非法输入后参数应恢复原值。');
assert(isfinite(str2double(get(ed,'String'))),'非法输入后编辑框应显示有限数值。');
set(ed,'String','1e9');   % 超范围，应被夹到上界
invoke(ed);
clamped=getappdata(f,'Parameters');
assert(clamped.Cd<=0.8+1e-12 && clamped.Cd>=0.08-1e-12,'超范围输入应被夹到允许区间。');
report.invalidInputRestored=true;
report.outOfRangeClamped=true;
assert(get(control('Param_windY'),'Min')==-15&&get(control('Param_windY'),'Max')==15);
assert(get(control('Param_density'),'Min')==.9&&get(control('Param_density'),'Max')==1.4);
report.expandedControls=true;

% ---------- 播放保持 ----------
tm=timerfind('Name','PenaltyRealtimePlayback');
assert(numel(tm)==1,'必须存在唯一的 PenaltyRealtimePlayback timer。');
play=control('PlayButton');
invoke(play);
assert(strcmp(tm.Running,'on'),'点击播放后 timer 应运行。');
ed=control('Edit_windY');
set(ed,'String','-2.00');
invoke(ed);
assert(strcmp(tm.Running,'on'),'播放中修改参数后应继续播放。');
sim=getappdata(f,'Scenario');
assert(getappdata(f,'CurrentTime')<sim.tEnd,'继续播放时不应直接跳到终点。');
assert(abs(getappdata(f,'Parameters').windY+2)<1e-12,'播放中改参未生效。');
invoke(play);   % 暂停
assert(strcmp(tm.Running,'off'),'再次点击应暂停 timer。');
report.playbackPreservedAfterEdit=true;

% ---------- 不合格候选：红色标示且不伪造结果 ----------
% 只用引擎返回的真实结果判定；若这些候选仍可行，则如实记录未触发，不伪造。
candidateWinds=[8,-8];
infeasibleFound=false;
for c=1:numel(candidateWinds)
    set(control('Edit_windY'),'String',sprintf('%.2f',candidateWinds(c)));
    invoke(control('Edit_windY'));
    sNow=getappdata(f,'Scenario');
    if isfield(sNow.fit,'feasible') && ~sNow.fit.feasible
        infeasibleFound=true;
        break;
    end
end
if infeasibleFound
    sNow=getappdata(f,'Scenario');
    assert(~sNow.fit.feasible,'不合格候选的 feasible 必须为 false。');
    cn=penalty_scenario('curve',sNow);
    assert(norm(cn.state(end,1:3)'-sNow.endpoint(1:3))<1e-6, ...
        '不合格候选也必须给出真实轨迹与末端，不得隐藏或伪造。');
    txt=get(st,'String');
    assert(~isempty(strfind(txt,'未满足观测容差')), ...
        '状态栏必须显著显示“该候选未满足观测容差”。');
    played=findobj(ax(1),'Tag','PlayedLine');
    fut=findobj(ax(1),'Tag','FutureLine');
    assert(all(abs(get(played,'Color')-[0.82 0.10 0.10])<1e-6), ...
        '不合格候选的已飞行轨迹必须为红色。');
    assert(all(abs(get(fut,'Color')-[0.82 0.10 0.10])<1e-6), ...
        '不合格候选的未来轨迹必须为红色。');
    fams=findall(ax(1),'Tag','FamilyLine');
    assert(~isempty(fams),'族线对象应存在（合格族可为空，但不加入当前候选）。');
    report.infeasibleMarkedRed=true;
else
    report.infeasibleMarkedRed=false;
    report.infeasibleNote='候选参数下引擎均返回可行解，未触发不合格分支（未伪造）。';
end

% ---------- 重播 ----------
invoke(control('ResetButton'));
resetSim=getappdata(f,'Scenario');
assert(abs(resetSim.params.tolerance-0.005)<1e-12,'复位后容差应回到默认值。');
assert(getappdata(f,'CurrentTime')==0,'复位后应停在 t=0。');
set(control('TimeSlider'),'Value',0.4*get(control('TimeSlider'),'Max'));
invoke(control('TimeSlider'));
assert(getappdata(f,'CurrentTime')>0,'重播前应先离开 0。');
invoke(control('ReplayButton'));
assert(getappdata(f,'CurrentTime')==0,'重播必须回到 t=0。');
assert(norm(getappdata(f,'Scenario').theta-resetSim.theta)<1e-12, ...
    '重播不应改变全程反演结果。');
report.replayReturnsToZero=true;

% ---------- 关闭清理 timer ----------
close(f);
assert(isempty(timerfind('Name','PenaltyRealtimePlayback')), ...
    '关闭界面后必须清理 PenaltyRealtimePlayback timer。');
report.timerClearedOnClose=true;

% ---------- 写回报告 ----------
report.passed=true;
outDir=fullfile(fileparts(mfilename('fullpath')),'data');
if ~exist(outDir,'dir'), mkdir(outDir); end
path=fullfile(outDir,'实时界面验证.json');
fid=fopen(path,'w','n','UTF-8');
assert(fid>=0,'无法写入验证报告。');
fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));
fclose(fid);
disp(report);

    function h=control(tag)
        h=findobj(f,'Tag',tag);
        assert(isscalar(h),['缺少控件 ' tag]);
    end
    function invoke(h)
        feval(h.Callback,h,[]);
    end
end

function closeIfLive(f)
if isgraphics(f),close(f);end
end
