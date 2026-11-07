function report=scenario_ui_checks()
%SCENARIO_UI_CHECKS 当前界面：一条完整计算线、无交点标签、改参后重求。
penalty_view;f=gcf;guard=onCleanup(@()closeIfLive(f));set(f,'Visible','off');
assert(numel(findall(f,'Type','axes'))==1);
assert(isempty(findall(f,'Tag','GoalIntersection'))&&isempty(findall(f,'Tag','GoalIntersectionLabel')));
assert(isempty(findall(f,'Tag','FutureLine')));
assert(control('FamilyToggle').Value==0);
line=control('PlayedLine');s=getappdata(f,'Scenario');curve=penalty_scenario('curve',s);
assert(norm(line.XData(:)-curve.state(:,1))<1e-12);
slider=control('TimeSlider');set(slider,'Value',.2);invoke(slider);
assert(norm(line.XData(:)-curve.state(:,1))<1e-12,'播放进度不改变整条计算线。');
old=s.theta;ed=control('Edit_windY');set(ed,'String','3.5');invoke(ed);
s=getappdata(f,'Scenario');assert(norm(s.theta(1:3)-old(1:3))>1e-5&&getappdata(f,'CurrentTime')==0);
assert(s.fit.feasible);
ed=control('Edit_Cd');before=s.params.Cd;set(ed,'String','NaN');invoke(ed);assert(getappdata(f,'Parameters').Cd==before);
play=control('PlayButton');invoke(play);timerH=timerfind('Name','PenaltyRealtimePlayback');assert(strcmp(timerH.Running,'on'));
set(control('Edit_windY'),'String','-8');invoke(control('Edit_windY'));assert(strcmp(timerH.Running,'on'));invoke(play);
set(control('Edit_Cd'),'String','.8');invoke(control('Edit_Cd'));s=getappdata(f,'Scenario');assert(~s.fit.feasible&&all(abs(line.Color-[.82 .1 .1])<1e-12));
invoke(control('ResetButton'));assert(control('FamilyToggle').Value==0);
assert(getappdata(f,'CurrentTime')==0);
close(f);assert(isempty(timerfind('Name','PenaltyRealtimePlayback')));
report=struct('passed',true,'singleCompleteCurve',true,'noPredictionAnnotation',true,'wholeCurveAtEveryTime',true,'parameterRefit',true,'playbackAndCleanup',true);
p=fullfile(fileparts(mfilename('fullpath')),'data','实时界面验证.json');fid=fopen(p,'w','n','UTF-8');fprintf(fid,'%s',jsonencode(report,'PrettyPrint',true));fclose(fid);disp(report);
 function h=control(tag),h=findall(f,'Tag',tag);assert(isscalar(h));end
 function invoke(h),feval(h.Callback,h,[]);end
end
function closeIfLive(f),if isgraphics(f),close(f);end,end
