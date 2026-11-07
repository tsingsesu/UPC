function report=scenario_ui_checks()
%SCENARIO_UI_CHECKS 实际回调集成检验：实时改参、模式、输入与资源生命周期。
penalty_view; f=gcf; set(f,'Visible','off');
cleanup=onCleanup(@()closeIfLive(f)); %#ok<NASGU>
assert(numel(findall(f,'Type','axes'))==1);
assert(strcmp(get(control('Param_speedScale'),'Enable'),'off'));
assert(strcmp(get(control('Param_mass'),'Enable'),'off'));
assert(~isempty(getappdata(f,'ParamListeners')));

sl=control('TimeSlider');set(sl,'Value',.2);invoke(sl);
original=getappdata(f,'Scenario');before=penalty_scenario('sample',original,.2);
speed=control('SpeedSelect');set(speed,'Value',4);invoke(speed);
tm=timerfind('Name','PenaltyRealtimePlayback');assert(numel(tm)==1);
set(tm,'TasksToExecute',200);
play=control('PlayButton');invoke(play);assert(strcmp(tm.Running,'on'));
wind=control('Param_windY');set(wind,'Value',12);tic;invoke(wind);report.liveEditLatency_s=toc;
changed=getappdata(f,'Scenario');now=penalty_scenario('sample',changed,.2);
assert(norm(now.state-before.state)<1e-9);
assert(abs(changed.params.windY-12)<1e-12);
assert(strcmp(tm.Running,'on'),'实时改参不应自动终止播放。');
assert(norm(changed.endpoint-original.endpoint)>.01);
invoke(play);assert(strcmp(tm.Running,'off'));
for t=[.05,.1,.15,.19]
    a=penalty_scenario('sample',original,t);b=penalty_scenario('sample',changed,t);
    assert(norm(a.state-b.state)<1e-9);
end

old=getappdata(f,'Parameters');edit=control('Edit_CdScale');set(edit,'String','NaN');invoke(edit);
new=getappdata(f,'Parameters');assert(old.CdScale==new.CdScale);
assert(isfinite(str2double(edit.String)));

mode=control('ModeSelect');set(mode,'Value',2);invoke(mode);
assert(strcmp(control('Param_speedScale').Enable,'on'));
mass=control('Param_mass');set(mass,'Value',.41);invoke(mass);
scale=control('Param_speedScale');set(scale,'Value',.8);invoke(scale);
free=getappdata(f,'Scenario');assert(strcmp(free.params.mode,'free'));
assert(abs(free.params.mass-.41)<1e-12&&abs(free.params.speedScale-.8)<1e-12);
assert(getappdata(f,'CurrentTime')==0&&strcmp(tm.Running,'off'));

set(mode,'Value',1);invoke(mode);fixed=getappdata(f,'Scenario');
assert(strcmp(fixed.params.mode,'conditioned'));
assert(abs(fixed.params.mass-.43)<1e-12&&fixed.params.speedScale==1);
for k=1:2
    q=penalty_scenario('sample',fixed,fixed.base.observations.times(k));
    assert(norm(q.state(1:3)-fixed.base.observations.positions(k,:)')<1e-7);
end

invoke(control('ResetButton'));reset=getappdata(f,'Scenario');
assert(norm(reset.endpoint-reset.reference.endpoint)<1e-7);
invoke(play);wait(tm);
assert(abs(getappdata(f,'CurrentTime')-getappdata(f,'Scenario').tEnd)<1e-9);
invoke(control('ReplayButton'));assert(getappdata(f,'CurrentTime')==0);
invoke(play);delete(f);assert(isempty(timerfind('Name','PenaltyRealtimePlayback')));
report.passed=true;report.singleAxes=true;report.causalLiveEdit=true;
report.modeResetAndObservationConstraint=true;report.invalidInputRejected=true;
path=fullfile(fileparts(mfilename('fullpath')),'data','实时界面验证.json');
fid=fopen(path,'w','n','UTF-8');assert(fid>=0);fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));fclose(fid);
disp(report);

    function h=control(tag)
        h=findobj(f,'Tag',tag);assert(isscalar(h),['缺少控件 ' tag]);
    end
    function invoke(h)
        feval(h.Callback,h,[]);
    end
end

function closeIfLive(f)
if isgraphics(f),close(f);end
end
