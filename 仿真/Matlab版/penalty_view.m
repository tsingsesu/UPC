function penalty_view(mode)
%PENALTY_VIEW 点球全程约束反演界面
%   penalty_view                交互界面（唯一可旋转/缩放三维轴 + 四个参数）
%   penalty_view('export')      导出 80 帧 GIF 与全景 PNG（默认全程参数）
%   penalty_view('demo')        从正式约束族选两组差异较大的合格成员演示
%   penalty_view('selftest')    轻量接口自检（不建界面）
%
%   参数（风、空气密度、阻力系数、观测容差）从出脚 t=0 起全程生效。
%   任何参数修改都重新做全程联合反演，不复用、不拼接旧历史轨迹。
%   动力学与拟合全部由 penalty_scenario 的稳定接口提供，本文件不实现第二套物理。

if nargin < 1 || isempty(mode), mode = 'interactive'; end
mode = lower(char(mode));
assert(any(strcmp(mode,{'interactive','export','demo','selftest'})), ...
    'mode 必须是 interactive / export / demo / selftest。');

thisFile = mfilename('fullpath'); outDir = fileparts(thisFile);
outDataDir = fullfile(outDir,'data'); outImgDir = fullfile(outDir,'img');

% ---------- 四个可调参数的范围 / 显示精度（API 一律用米、m/s、kg/m^3） ----------
PRM = struct( ...
    'names', {{'windY','density','Cd','tolerance'}}, ...
    'label', {{'横风 windY（m/s）','空气密度（kg/m^3）','阻力系数 Cd','观测容差（m）'}}, ...
    'lo',    [-15, .90, .08, 0.001], ...
    'hi',    [15, 1.40, .80, 0.020], ...
    'dp',    [ 2, 3, 6, 4], ...
    'step',  {{[0.002 0.05],[0.0005 0.01],[0.0005 0.01],[0.0001 0.01]}});

% ---------- 引擎 ----------
p = withDefaults(penalty_scenario('defaults'));
sim = penalty_scenario('start',p);
assertSim(sim);

if strcmp(mode,'selftest')
    runSelftest();
    return;
end

if any(strcmp(mode,{'export','demo'}))
    if ~exist(outImgDir,'dir'), mkdir(outImgDir); end
    if ~exist(outDataDir,'dir'), mkdir(outDataDir); end
end

% ---------- 运行时状态 ----------
playing = false; PLAY_SPEED = 0.2; t0 = 0; tClock = []; timerH = [];
curveT = []; curveS = [];
family = []; familyLoaded = false; familyTol = NaN;
familyQualified = 0; familyDrawn = 0; familyError = '';
statusMsg = ''; statusWarn = false;
fnUI = 'Microsoft YaHei';
scaleF = 0.25;             % 力箭头：0.25 m / N

% 球门几何（与正式基准一致）
xPlane = 11; goalHW = 3.66; goalH = 2.44;

% ---------- 图窗与唯一三维轴 ----------
fig = figure('Name','点球全程约束反演','Color','w','Units','pixels', ...
    'Position',[40 40 1440 860],'NumberTitle','off');
if any(strcmp(mode,{'export','demo'}))
    set(fig,'Visible','off');
end
ax3 = axes(fig,'Units','normalized','Position',[0.395 0.205 0.575 0.720], ...
    'Tag','Trajectory3D');
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on'); view(ax3,-32,18);
set(ax3,'FontName',fnUI);

% 球门（与真实尺寸一致）
patch(ax3,[9 11.4 11.4 9],[-0.4 -0.4 3.8 3.8],[0 0 0 0],[0.42 0.72 0.42], ...
    'FaceColor',[0.42 0.72 0.42],'EdgeColor','none','FaceAlpha',0.10, ...
    'HandleVisibility','off');
plot3(ax3,xPlane*[1 1], goalHW*[1 1],[0 goalH],'-','Color',[0.12 0.12 0.12], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],-goalHW*[1 1],[0 goalH],'-','Color',[0.12 0.12 0.12], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],[-goalHW goalHW],goalH*[1 1],'-','Color',[0.12 0.12 0.12], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],[-goalHW goalHW],[0 0],'-','Color',[0.12 0.12 0.12], ...
    'LineWidth',1.0,'HandleVisibility','off');

% 半透明门面表示x=11的门口区域，交点仍由求解器事件定位。
patch(ax3,11*ones(1,4),[-goalHW goalHW goalHW -goalHW],[0 0 goalH goalH], ...
    [.20 .65 .85],'FaceAlpha',.045,'EdgeColor','none','HandleVisibility','off','Tag','GoalPlane');
hFamilyGoal=plot3(ax3,nan,nan,nan,'.','Color',[.40 .40 .40],'MarkerSize',9, ...
    'HandleVisibility','off','Tag','FamilyGoalIntersections');
hGoalPoint=plot3(ax3,nan,nan,nan,'d','Color',[.55 .12 .70], ...
    'MarkerFaceColor','none','MarkerSize',18,'LineWidth',2, ...
    'DisplayName','当前预测交点','Tag','GoalIntersection');
hGoalCross=plot3(ax3,nan,nan,nan,'-','Color',[.55 .12 .70], ...
    'LineWidth',1.4,'HandleVisibility','off','Tag','GoalIntersectionCross');
hGoalText=text(ax3,0,0,0,'','FontName',fnUI,'FontSize',9, ...
    'Color',[.45 .08 .60],'BackgroundColor','w','Margin',3,'Clipping','on', ...
    'HandleVisibility','off','Tag','GoalIntersectionLabel');

% 灰色合格族（最多画 80 条，均匀抽样）
maxFamLines = 80;
famLines = gobjects(maxFamLines,1);
for k = 1:maxFamLines
    famLines(k) = plot3(ax3,nan,nan,nan,'-','Color',[0.72 0.72 0.72], ...
        'LineWidth',0.5,'HandleVisibility','off','Tag','FamilyLine');
end
hFamLegend = plot3(ax3,nan,nan,nan,'-','Color',[0.72 0.72 0.72], ...
    'LineWidth',0.8,'DisplayName','灰色族及其交点');

% 已飞行 / 未来轨迹
hPlayed = plot3(ax3,nan,nan,nan,'-','Color',[0 0.35 0.85], ...
    'LineWidth',2.2,'DisplayName','已飞行','Tag','PlayedLine');
hFuture = plot3(ax3,nan,nan,nan,'--','Color',[0.95 0.55 0.08], ...
    'LineWidth',1.6,'DisplayName','未来','Tag','FutureLine');

% 题给观测 + 当前拟合位置标记
hObs = plot3(ax3,nan,nan,nan,'ko','MarkerFaceColor','k','MarkerSize',6, ...
    'DisplayName','题给观测','Tag','Observed');
hFitObs = plot3(ax3,nan,nan,nan,'ks','MarkerSize',8,'LineWidth',1.2, ...
    'DisplayName','拟合观测位置','Tag','FitObserved');
hCur = plot3(ax3,nan,nan,nan,'kd','MarkerFaceColor','k','MarkerSize',5, ...
    'LineWidth',1.0,'DisplayName','当前拟合位置','Tag','CurrentFit');

% 足球（固定真实半径 .11 m）
[Xb,Yb,Zb] = sphere(20);
hBall = surf(ax3,Xb,Yb,Zb,'FaceColor',[0.97 0.97 0.97], ...
    'EdgeColor',[0.25 0.25 0.25],'LineWidth',0.3,'HandleVisibility','off');

% 三个力箭头（红重力 / 绿阻力 / 蓝 Magnus），初始隐藏
hqG = quiver3(ax3,0,0,0,0,0,0,'Color',[0.85 0.2 0.2],'LineWidth',1.5, ...
    'MaxHeadSize',0.8,'HandleVisibility','off');
hqD = quiver3(ax3,0,0,0,0,0,0,'Color',[0.2 0.65 0.2],'LineWidth',1.5, ...
    'MaxHeadSize',0.8,'HandleVisibility','off');
hqM = quiver3(ax3,0,0,0,0,0,0,'Color',[0.2 0.3 0.9],'LineWidth',1.5, ...
    'MaxHeadSize',0.8,'HandleVisibility','off');
set([hqG hqD hqM],'Visible','off');

% demo 对比专用（平时隐藏）
hCmp = gobjects(1,2);
hCmpTxt = gobjects(1,2);
hCmpGoal=gobjects(1,2);
hCmp(1) = plot3(ax3,nan,nan,nan,'-','Color',[0 0.35 0.85],'LineWidth',2.0, ...
    'HandleVisibility','off');
hCmp(2) = plot3(ax3,nan,nan,nan,'-','Color',[0.92 0.45 0.05],'LineWidth',2.0, ...
    'HandleVisibility','off');
hCmpTxt(1) = text(ax3,0,0,0,'第 1 组全程参数','FontName',fnUI,'FontSize',9, ...
    'Color',[0 0.35 0.85],'FontWeight','bold','Visible','off');
hCmpTxt(2) = text(ax3,0,0,0,'第 2 组全程参数','FontName',fnUI,'FontSize',9, ...
    'Color',[0.92 0.45 0.05],'FontWeight','bold','Visible','off');
set(hCmp,'Visible','off');
for kk=1:2
 hCmpGoal(kk)=plot3(ax3,nan,nan,nan,'d','MarkerSize',9,'LineWidth',1.5, ...
  'Color',hCmp(kk).Color,'MarkerFaceColor',hCmp(kk).Color, ...
  'Visible','off','HandleVisibility','off','Tag','ComparisonGoalIntersection');
end

xlabel(ax3,'x（m，指向球门）'); ylabel(ax3,'y（m，横向）'); zlabel(ax3,'z（m，高度）');
title(ax3,'点球全程约束反演（灰=合格参数族；蓝=已飞行；橙=未来）', ...
    'FontName',fnUI,'FontSize',11);
lg = legend(ax3,[hFamLegend hObs hFitObs hPlayed hFuture hCur hGoalPoint], ...
    'Location','northwest','FontSize',8,'Box','off');
set(lg,'AutoUpdate','off');

% ---------- 左侧参数面板 ----------
uicontrol(fig,'Style','text','Units','normalized','Position',[0.015 0.947 0.30 0.038], ...
    'String','全程约束反演参数（改动后从 t=0 重新反演）','HorizontalAlignment','left', ...
    'FontName',fnUI,'FontSize',10,'FontWeight','bold');

hS = struct(); hE = struct();
for ii = 1:numel(PRM.names)
    nm = PRM.names{ii};
    y = 0.895 - (ii-1)*0.056;
    uicontrol(fig,'Style','text','Units','normalized','Position',[0.015 y 0.105 0.042], ...
        'String',PRM.label{ii},'HorizontalAlignment','left','FontName',fnUI,'FontSize',9);
    hS.(nm) = uicontrol(fig,'Style','slider','Units','normalized', ...
        'Position',[0.125 y+0.008 0.110 0.030],'Min',PRM.lo(ii),'Max',PRM.hi(ii), ...
        'Value',p.(nm),'Tag',['Param_' nm],'SliderStep',PRM.step{ii}, ...
        'Callback',@(~,~)cbSlider(nm), ...
        'ButtonDownFcn',@(~,~)setStatus('松开滑块后重新反演。'));
    hE.(nm) = uicontrol(fig,'Style','edit','Units','normalized', ...
        'Position',[0.242 y+0.004 0.073 0.038],'String',fmt(p.(nm),PRM.dp(ii)), ...
        'Tag',['Edit_' nm],'Callback',@(~,~)cbEdit(nm), ...
        'FontName','Consolas','FontSize',9);
end

uicontrol(fig,'Style','text','Units','normalized','Position',[0.015 0.665 0.30 0.030], ...
    'String','—— 全程联合反演结果 ——','HorizontalAlignment','left', ...
    'FontName',fnUI,'FontSize',9,'FontWeight','bold');

hInfo = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.015 0.205 0.305 0.455],'String','', ...
    'HorizontalAlignment','left', ...
    'FontName',fnUI,'FontSize',9,'Tag','InfoText');

hFamNote = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.015 0.145 0.305 0.055],'String','', ...
    'HorizontalAlignment','left', ...
    'FontName',fnUI,'FontSize',8,'Tag','FamilyNote');

% ---------- 底部控制 ----------
hPlay = uicontrol(fig,'Style','pushbutton','Units','normalized', ...
    'Position',[0.335 0.095 0.055 0.040],'String','播放','Tag','PlayButton', ...
    'Callback',@(~,~)cbPlay(),'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','pushbutton','Units','normalized', ...
    'Position',[0.393 0.095 0.055 0.040],'String','重播','Tag','ReplayButton', ...
    'Callback',@(~,~)cbReplay(),'FontName',fnUI,'FontSize',9);
hTime = uicontrol(fig,'Style','slider','Units','normalized', ...
    'Position',[0.451 0.100 0.200 0.030],'Min',0,'Max',max(sim.tEnd,1e-6), ...
    'Value',0,'SliderStep',[0.002 0.05],'Tag','TimeSlider', ...
    'Callback',@(~,~)cbScrub(), ...
    'ButtonDownFcn',@(~,~)setStatus('松开滑块后定位时间。'));
hView = uicontrol(fig,'Style','popupmenu','Units','normalized', ...
    'Position',[0.655 0.095 0.105 0.040],'String',{'全景','末段放大'},'Value',1, ...
    'Tag','ViewSelect','Callback',@(~,~)cbView(),'FontName',fnUI,'FontSize',9);
hFamilyToggle = uicontrol(fig,'Style','checkbox','Units','normalized', ...
    'Position',[0.765 0.097 0.085 0.036],'String','族显示','Value',1, ...
    'Tag','FamilyToggle','Callback',@(~,~)cbFamily(),'FontName',fnUI,'FontSize',9);
hForceToggle = uicontrol(fig,'Style','checkbox','Units','normalized', ...
    'Position',[0.855 0.097 0.065 0.036],'String','力箭头','Value',0, ...
    'Tag','ForceToggle','Callback',@(~,~)cbForce(),'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','pushbutton','Units','normalized', ...
    'Position',[0.925 0.095 0.062 0.040],'String','复位','Tag','ResetButton', ...
    'Callback',@(~,~)cbReset(),'FontName',fnUI,'FontSize',9);

hTimeReadout = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.335 0.058 0.305 0.030],'String','', ...
    'HorizontalAlignment','left','FontName',fnUI,'FontSize',8,'Tag','TimeReadout');
uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.015 0.052 0.305 0.036],'String','松开后重新反演；宽范围为探索情景', ...
    'HorizontalAlignment','left','FontName',fnUI,'FontSize',8, ...
    'ForegroundColor',[0.35 0.35 0.35],'Tag','SliderHint');

hStatus = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.015 0.005 0.972 0.045],'String','', ...
    'HorizontalAlignment','left','FontName',fnUI,'FontSize',10, ...
    'Tag','StatusText');

% ---------- appdata ----------
setappdata(fig,'Scenario',sim);
setappdata(fig,'Parameters',p);
setappdata(fig,'CurrentTime',0);
setappdata(fig,'Busy',false);
setappdata(fig,'PlaybackTimer',[]);

set(fig,'CloseRequestFcn',@(~,~)cbClose());
set(fig,'DeleteFcn',@(~,~)onDelete());

% ---------- 初始化显示 ----------
refreshCurve();
refreshObservations();
refreshFamily(true);
syncControls();
resetTimeZero();
applyStyles();
applyView();
drawAt(0);
if sim.fit.feasible
    setStatus(sprintf('就绪：全程反演 T = %.3f s，最大残差 %.2f mm ≤ 容差 %.2f mm。', ...
        sim.tEnd, sim.fit.maxResidual*1000, sim.params.tolerance*1000));
else
    setStatus(candidateStatus(sim),true);
end

if strcmp(mode,'export')
    runExport();
    return;
elseif strcmp(mode,'demo')
    runDemo();
    return;
end

% ---------- 交互播放 ----------
timerH = timer('ExecutionMode','fixedSpacing','Period',0.04,'BusyMode','drop', ...
    'TimerFcn',@(~,~)cbTick(),'Name','PenaltyRealtimePlayback');
setappdata(fig,'PlaybackTimer',timerH);
rotate3d(fig,'on');
fprintf('[penalty_view] 就绪：T=%.6f s，状态 %s，最大残差 %.6f m。\n', ...
    sim.tEnd, sim.status, sim.fit.maxResidual);

% ================= 界面回调 =================
    function cbSlider(name)
        if busyNow(), return; end
        applyParam(name,get(hS.(name),'Value'));
    end

    function cbEdit(name)
        if busyNow(), return; end
        raw = get(hE.(name),'String');
        val = str2double(raw);
        if isnan(val) || ~isfinite(val)
            set(hE.(name),'String',fmt(p.(name),PRM.dp(paramIndex(name))));
            setStatus(sprintf('%s 输入无效（非数字/NaN/Inf），已恢复原值。',PRM.label{paramIndex(name)}),true);
            return;
        end
        applyParam(name,val);
    end

    function applyParam(name,val)
        k = paramIndex(name);
        oldP = p; oldSim = sim; oldTime = getappdata(fig,'CurrentTime');
        wasPlaying = playing;
        stopPlayback();

        setappdata(fig,'Busy',true);
        busyCleanup = onCleanup(@()setappdata(fig,'Busy',false)); %#ok<NASGU>
        try
            assert(isnumeric(val) && isscalar(val) && isfinite(val),'参数必须是有限标量。');
            val = min(max(val,PRM.lo(k)),PRM.hi(k));
            val = round(val*10^PRM.dp(k))/10^PRM.dp(k);

            p.(name) = val;
            setappdata(fig,'Parameters',p);
            syncControls();

            % 全程重算：从 t=0 的联合拟合开始，绝不拼接旧历史。
            sim = penalty_scenario('start',p);
            assertSim(sim);
            setappdata(fig,'Scenario',sim);

            refreshCurve();
            refreshObservations();
            if strcmp(name,'tolerance'), refreshFamily(true); end
            resetTimeZero();
            applyStyles();
            drawAt(0);

            if wasPlaying
                startPlayback();
                setStatus(sprintf('%s 已更新：全程从 t=0 重新反演，并继续播放。',PRM.label{k}));
            else
                setStatus(sprintf('%s 已更新：全程从 t=0 重新反演，停在 t=0。',PRM.label{k}));
            end
            if ~sim.fit.feasible
                setStatus(candidateStatus(sim),true);
            end
        catch err
            p = oldP; sim = oldSim;
            setappdata(fig,'Scenario',sim);
            setappdata(fig,'Parameters',p);
            syncControls();
            refreshCurve(); refreshObservations();
            set(hTime,'Max',max(sim.tEnd,1e-6));
            drawAt(min(oldTime,sim.tEnd));
            applyStyles();
            if wasPlaying
                t0 = getappdata(fig,'CurrentTime'); tClock = tic; playing = true;
                if ~isempty(timerH) && isvalid(timerH), start(timerH); end
                if isgraphics(hPlay), set(hPlay,'String','暂停'); end
            end
            setStatus(['反演失败，已恢复原参数：' err.message],true);
        end
    end

    function cbPlay()
        if busyNow(), return; end
        if playing
            stopPlayback();
            setStatus(sprintf('已暂停于 t = %.3f s。',getappdata(fig,'CurrentTime')));
            return;
        end
        if getappdata(fig,'CurrentTime') >= sim.tEnd - 1e-9
            drawAt(0);
        end
        startPlayback();
        setStatus('播放中：环境参数从出脚起全程生效。');
        if ~sim.fit.feasible, setStatus(candidateStatus(sim),true); end
    end

    function startPlayback()
        if sim.tEnd <= 0, return; end
        t0 = getappdata(fig,'CurrentTime');
        tClock = tic; playing = true;
        if ~isempty(timerH) && isvalid(timerH), start(timerH); end
        if isgraphics(hPlay), set(hPlay,'String','暂停'); end
    end

    function stopPlayback()
        playing = false;
        stopOnly(timerH);
        if isgraphics(hPlay), set(hPlay,'String','播放'); end
    end

    function cbTick()
        if ~playing || busyNow(), return; end
        if isempty(tClock), tClock = tic; end
        tn = t0 + PLAY_SPEED*toc(tClock);
        if tn >= sim.tEnd
            drawAt(sim.tEnd);
            playing = false;
            stopOnly(timerH);
            if isgraphics(hPlay), set(hPlay,'String','播放'); end
            setStatus(sprintf('播放结束：t = %.3f s，末端状态 %s。',sim.tEnd,statusChinese(sim.status)));
            if ~sim.fit.feasible, setStatus(candidateStatus(sim),true); end
            return;
        end
        drawAt(tn);
    end

    function cbReplay()
        if busyNow(), return; end
        stopPlayback();
        drawAt(0);
        setStatus('已回到 t=0，点击播放重看当前参数下的射门。');
    end

    function cbScrub()
        if busyNow(), return; end
        stopPlayback();
        tn = get(hTime,'Value');
        drawAt(tn);
        setStatus(sprintf('时间已定位到 t = %.3f s（拖动不改变全程反演结果）。',tn));
    end

    function cbView()
        applyView();
        setStatus('视角已切换。');
    end

    function cbFamily()
        updateFamilyVisibility();
        updateFamilyNote();
        if isempty(family)
            setStatus('正式约束族数据尚未生成，灰色族为空。');
        else
            setStatus(sprintf('灰色族显示已切换：合格 %d 条，实际绘制 %d 条（不是概率）。', ...
                familyQualified, familyDrawn));
        end
    end

    function cbForce()
        updateForceVisibility();
        if get(hForceToggle,'Value')
            setStatus('力箭头已显示（红=重力，绿=阻力，蓝=Magnus；0.25 m/N）。');
        else
            setStatus('力箭头已隐藏。');
        end
    end

    function cbReset()
        if busyNow(), return; end
        oldP = p; oldSim = sim;
        stopPlayback();
        setappdata(fig,'Busy',true);
        busyCleanup = onCleanup(@()setappdata(fig,'Busy',false)); %#ok<NASGU>
        try
            p = withDefaults(penalty_scenario('defaults'));
            sim = penalty_scenario('start',p);
            assertSim(sim);
            setappdata(fig,'Scenario',sim);
            setappdata(fig,'Parameters',p);
            set(hView,'Value',1);
            set(hForceToggle,'Value',0);
            set(hFamilyToggle,'Value',1);
            refreshCurve(); refreshObservations(); refreshFamily(true);
            syncControls(); resetTimeZero();
            applyView(); updateForceVisibility(); updateFamilyVisibility();
            applyStyles(); drawAt(0);
            setStatus('已复位：默认参数、t=0、暂停；全程反演已重新完成。');
        catch err
            p = oldP; sim = oldSim;
            setappdata(fig,'Scenario',sim);
            setappdata(fig,'Parameters',p);
            syncControls();
            setStatus(['复位失败，已保留原状态：' err.message],true);
        end
    end

    function cbClose()
        stopPlayback();
        stopTimer(timerH);
        if isgraphics(fig), delete(fig); end
    end

    function onDelete()
        stopPlayback();
        stopTimer(timerH);
    end

    % ================= 显示 =================
    function refreshCurve()
        c = penalty_scenario('curve',sim);
        t = c.t(:); s = c.state;
        if size(s,1) ~= numel(t), s = s.'; end
        assert(size(s,2) >= 6,'curve 状态必须为 N×6。');
        curveT = t; curveS = s;
        refreshGoalIntersection();
    end

    function refreshGoalIntersection()
        % 只标出真正到达x=11的事件；先触地的端点不是球门交点。
        ep=sim.endpoint(:);
        crossing=strcmp(sim.status,'plane')&&numel(ep)>=3&&all(isfinite(ep(1:3)))&&abs(ep(1)-xPlane)<1e-7;
        if ~crossing
            set([hGoalPoint hGoalCross hGoalText],'Visible','off');
            set(hGoalPoint,'XData',nan,'YData',nan,'ZData',nan);
            return;
        end
        color=[.55 .12 .70];note='预测交点';
        if ~sim.fit.feasible,color=[.82 .10 .10];note='未通过筛选的候选交点';end
        set(hGoalPoint,'XData',ep(1),'YData',ep(2),'ZData',ep(3),'Color',color,'MarkerFaceColor','none','Visible','on');
        a=.23;
        set(hGoalCross,'XData',ep(1)*ones(1,5), ...
            'YData',[ep(2)-a ep(2)+a nan ep(2) ep(2)], ...
            'ZData',[ep(3) ep(3) nan ep(3)-a ep(3)+a],'Color',color,'Visible','on');
        set(hGoalText,'Position',[ep(1)-.65 ep(2)+.15 ep(3)+.26], ...
            'String',sprintf('%s\n(11, %.3f, %.3f) m',note,ep(2),ep(3)), ...
            'Color',color,'Visible','on');
    end

    function refreshObservations()
        obsT = sim.base.observations.times(:);
        obs = sim.base.observations.positions;
        if size(obs,1) ~= numel(obsT), obs = obs.'; end
        fitPos = nan(numel(obsT),3);
        for kk = 1:numel(obsT)
            try
                q = penalty_scenario('sample',sim,obsT(kk));
                fitPos(kk,:) = q.state(1:3)';
            catch
            end
        end
        set(hObs,'XData',obs(:,1),'YData',obs(:,2),'ZData',obs(:,3));
        set(hFitObs,'XData',fitPos(:,1),'YData',fitPos(:,2),'ZData',fitPos(:,3));
    end

    function drawAt(tNow)
        tNow = max(0,min(tNow,sim.tEnd));
        sm = penalty_scenario('sample',sim,tNow);
        s = sm.state(:);
        S = s(1:3)';

        rad = 0.11;
        if isfield(sim.params,'radius') && isfinite(sim.params.radius) && sim.params.radius > 0
            rad = sim.params.radius;
        end
        set(hBall,'XData',S(1)+rad*Xb,'YData',S(2)+rad*Yb,'ZData',S(3)+rad*Zb);

        k = find(curveT <= tNow+1e-12);
        if isempty(k)
            P = [S;S];
        else
            P = curveS(k,1:3);
            if curveT(k(end)) < tNow-1e-12, P = [P;S]; end
        end
        kf = find(curveT > tNow+1e-12);
        Ftr = [S;curveS(kf,1:3)];
        set(hPlayed,'XData',P(:,1),'YData',P(:,2),'ZData',P(:,3));
        set(hFuture,'XData',Ftr(:,1),'YData',Ftr(:,2),'ZData',Ftr(:,3));
        set(hCur,'XData',S(1),'YData',S(2),'ZData',S(3));

        F = sm.forces;
        if ~isfield(F,'gravity'), F.gravity = [0;0;0]; end
        if ~isfield(F,'drag'),    F.drag    = [0;0;0]; end
        if ~isfield(F,'magnus'),  F.magnus  = [0;0;0]; end
        set(hqG,'XData',S(1),'YData',S(2),'ZData',S(3), ...
            'UData',scaleF*F.gravity(1),'VData',scaleF*F.gravity(2),'WData',scaleF*F.gravity(3));
        set(hqD,'XData',S(1),'YData',S(2),'ZData',S(3), ...
            'UData',scaleF*F.drag(1),'VData',scaleF*F.drag(2),'WData',scaleF*F.drag(3));
        set(hqM,'XData',S(1),'YData',S(2),'ZData',S(3), ...
            'UData',scaleF*F.magnus(1),'VData',scaleF*F.magnus(2),'WData',scaleF*F.magnus(3));

        set(hTime,'Value',tNow);
        setappdata(fig,'CurrentTime',tNow);
        set(hTimeReadout,'String',sprintf('t = %.3f s  ／  总时长 %.3f s',tNow,sim.tEnd));
        renderInfo(tNow,sm);
    end

    function renderInfo(tNow,sm)
        s = sm.state(:);
        th = sim.theta(:);
        if numel(th) < 6, th = zeros(6,1); end
        v0 = th(1:3);
        sp = NaN;
        if isfield(sim,'initialSpeed') && isfinite(sim.initialSpeed), sp = sim.initialSpeed; end
        if ~isfinite(sp), sp = norm(v0); end
        E = NaN;
        if isfield(sim.fit,'kineticEnergy') && isfinite(sim.fit.kineticEnergy)
            E = sim.fit.kineticEnergy;
        end
        if ~isfinite(E), E = 0.5*0.43*sp^2; end
        tol = 0.005;
        if isfield(sim.params,'tolerance'), tol = sim.params.tolerance; end
        mr = NaN;
        if isfield(sim.fit,'maxResidual'), mr = sim.fit.maxResidual; end
        feas = isfield(sim.fit,'feasible') && sim.fit.feasible;
        if isfinite(mr)
            residTxt = sprintf('%.2f mm',mr*1000);
        else
            residTxt = '—';
        end
        residList = '—';
        if isfield(sim.fit,'residual') && ~isempty(sim.fit.residual)
            R = sim.fit.residual;
            rn = max(abs(R),[],2)*1000;
            parts = arrayfun(@(v)sprintf('%.2f',v),rn(:),'UniformOutput',false);
            residList = [strjoin(parts,' / ') ' mm'];
        end
        ep = sim.endpoint(:);
        if numel(ep) < 3, ep = [NaN;NaN;NaN]; end
        if feas
            feasTxt = '满足观测与物理筛选';
        else
            feasTxt = candidateStatus(sim);
        end
        if isfield(sim,'isOpening') && sim.isOpening
            openTxt = '进门';
        else
            openTxt = '未进门';
        end
        w = sm.wind(:);
        if numel(w) < 3, w = [0;0;0]; end
        CdNow = NaN;
        if isfield(sm,'Cd'), CdNow = sm.Cd; end
        lines = { ...
            sprintf('t = %.3f s ｜ 速度 %.2f m/s ｜ 位置 (%.2f, %.2f, %.2f) m', ...
                tNow,norm(s(4:6)),s(1),s(2),s(3)), ...
            '—— 起脚条件（全程联合反演） ——', ...
            sprintf('v0 = (%+.2f, %+.2f, %+.2f) m/s',v0(1),v0(2),v0(3)), ...
            sprintf('起脚速度 %.2f m/s ｜ 起脚动能 %.2f J（不是脚的作用力）',sp,E), ...
            sprintf('指定 Cd = %.4f ｜ by = %+.3f ｜ bz = %+.3f',th(4),th(5),th(6)), ...
            sprintf('最大侧力系数 %.3f / 上限 %.2f',sim.admissibility.maxTransverseCoefficient,sim.admissibility.limits.maxTransverseCoefficient), ...
            '—— 拟合质量与末端 ——', ...
            sprintf('最大残差 %s ｜ 容差 %.2f mm',residTxt,tol*1000), ...
            sprintf('各时刻最大坐标残差：%s',residList), ...
            sprintf('末端（%s）：(%+.2f, %+.2f, %+.2f) m',statusChinese(sim.status),ep(1),ep(2),ep(3)), ...
            sprintf('合格状态：%s ｜ %s',feasTxt,openTxt), ...
            sprintf('当前风 (%+.2f, %+.2f, %+.2f) m/s ｜ 有效 Cd %.4f',w(1),w(2),w(3),CdNow), ...
            '环境从出脚 t=0 起全程生效。', ...
            '灰色族为相容解的采样子集，不是置信区间。'};
        set(hInfo,'String',lines);
    end

    function applyStyles()
        feas = isfield(sim.fit,'feasible') && sim.fit.feasible;
        if feas
            set(hPlayed,'Color',[0 0.35 0.85]);
            set(hFuture,'Color',[0.95 0.55 0.08]);
            set(ax3.Title,'String','点球全程约束反演（灰=合格参数族；蓝=已飞行；橙=未来）', ...
                'Color',[0.1 0.1 0.1]);
        else
            set(hPlayed,'Color',[0.82 0.10 0.10]);
            set(hFuture,'Color',[0.82 0.10 0.10]);
            set(ax3.Title,'String',candidateStatus(sim), ...
                'Color',[0.75 0 0]);
        end
    end

    function applyView()
        if get(hView,'Value') == 1
            xr = [-0.8 11.8]; yr = [-4.3 4.3]; zr = [-0.05 3.3];
        else
            xr = [9.0 11.3]; yr = [1.9 3.1]; zr = [1.5 2.5];
        end
        xlim(ax3,xr); ylim(ax3,yr); zlim(ax3,zr);
        pbaspect(ax3,[diff(xr) diff(yr) diff(zr)]);
    end

    function updateForceVisibility()
        onoff = 'off';
        if get(hForceToggle,'Value'), onoff = 'on'; end
        set([hqG hqD hqM],'Visible',onoff);
    end

    function syncControls()
        for jj = 1:numel(PRM.names)
            nm = PRM.names{jj};
            set(hS.(nm),'Value',p.(nm));
            set(hE.(nm),'String',fmt(p.(nm),PRM.dp(jj)));
        end
    end

    function resetTimeZero()
        set(hTime,'Min',0,'Max',max(sim.tEnd,1e-6),'Value',0);
        setappdata(fig,'CurrentTime',0);
        set(hTimeReadout,'String',sprintf('t = 0.000 s  ／  总时长 %.3f s',sim.tEnd));
    end

    function setStatus(msg,warn)
        if nargin < 2, warn = false; end
        statusMsg = msg; statusWarn = warn;
        infeasible = isfield(sim,'fit') && isfield(sim.fit,'feasible') && ~sim.fit.feasible;
        if infeasible
            txt = ['【该候选未满足观测容差】 ' msg];
            col = [0.75 0 0]; fw = 'bold';
        elseif warn
            txt = msg; col = [0.75 0 0]; fw = 'bold';
        else
            txt = msg; col = [0.15 0.15 0.15]; fw = 'normal';
        end
        if isgraphics(hStatus)
            set(hStatus,'String',txt,'ForegroundColor',col,'FontWeight',fw);
        end
    end

    % ================= 约束族 =================
    function refreshFamily(force)
        if nargin < 1, force = false; end
        if force || ~familyLoaded || ~isfinite(familyTol) || abs(familyTol-p.tolerance) > 1e-12
            familyError = '';
            family = [];
            try
                fam = penalty_scenario('family',p.tolerance);
                if isstruct(fam), family = fam; end
            catch err
                familyError = err.message;
            end
            if ~isempty(family) && isfield(family,'accepted')
                family = family([family.accepted] ~= 0);
            end
            familyTol = p.tolerance;
            familyLoaded = true;
        end
        drawFamily();
    end

    function drawFamily()
        famMax = 0;
        if ~isempty(family), famMax = numel(family); end
        familyQualified = famMax;
        familyDrawn = 0;
        if famMax > 0
            idx = unique(round(linspace(1,famMax,min(80,famMax))));
            familyDrawn = numel(idx);
        else
            idx = [];
        end
        crossings=zeros(0,3);
        for kk = 1:maxFamLines
            if kk <= familyDrawn
                member=family(idx(kk));tr=member.trajectory;
                if strcmp(member.status,'plane')&&abs(member.endpoint(1)-xPlane)<1e-7
                    crossings(end+1,:)=member.endpoint(1:3)'; %#ok<AGROW>
                end
                if size(tr,2) >= 3
                    set(famLines(kk),'XData',tr(:,1),'YData',tr(:,2),'ZData',tr(:,3));
                else
                    set(famLines(kk),'XData',NaN,'YData',NaN,'ZData',NaN);
                end
            else
                set(famLines(kk),'XData',NaN,'YData',NaN,'ZData',NaN);
            end
        end
        set(hFamilyGoal,'XData',crossings(:,1),'YData',crossings(:,2),'ZData',crossings(:,3));
        updateFamilyVisibility();
        updateFamilyNote();
    end

    function updateFamilyVisibility()
        vis = 'off';
        if get(hFamilyToggle,'Value'), vis = 'on'; end
        set(hFamilyGoal,'Visible',vis);
        for kk = 1:familyDrawn
            set(famLines(kk),'Visible',vis);
        end
    end

    function updateFamilyNote()
        vis = get(hFamilyToggle,'Value') ~= 0;
        if isempty(family)
            if isempty(familyError)
                txt = '灰色族：正式约束族数据尚未生成，暂无灰色轨迹（不是置信区间/概率）。';
            else
                txt = ['灰色族不可用：' familyError];
            end
        else
            txt = sprintf('论文网格族：合格 %d 条，实际绘制 %d 条（等间隔抽样）；不是置信区间，也不是概率。', ...
                familyQualified,familyDrawn);
            if ~vis, txt = [txt '（已隐藏）']; end
        end
        set(hFamNote,'String',txt);
    end

    % ================= 导出 =================
    function runExport()
        gifPath = fullfile(outImgDir,'点球物理仿真.gif');
        pngPath = fullfile(outImgDir,'点球仿真全景.png');
        nF = 80; gifDur = 4.5;
        for kk = 1:nF
            tNow = sim.tEnd*(kk-1)/(nF-1);
            drawAt(tNow);
            setStatus(sprintf('导出全景：第 %d / %d 帧。',kk,nF));
            drawnow;
            writeGifFrame(fig,gifPath,kk,gifDur/nF);
        end
        drawAt(sim.tEnd);
        setStatus('导出全景：完成。');
        drawnow;
        fr = getframe(fig);
        imwrite(frame2im(fr),pngPath);

        curve=penalty_scenario('curve',sim);
        T=array2table([curve.t curve.state],'VariableNames',{'t_s','x_m','y_m','z_m','vx_mps','vy_mps','vz_mps'});
        writetable(T,fullfile(outDataDir,'点球仿真时序.csv'));
        J = buildMeta(sim,p,'export',struct('gif',gifPath,'png',pngPath));
        writeJson(fullfile(outDataDir,'实时参数导出.json'),J);
        close(fig);
        fprintf('[penalty_view] export 完成：%s，%s。\n',gifPath,pngPath);
    end

    function runDemo()
        if isempty(family)
            family = loadFamilyDirect(p.tolerance);
        end
        if isempty(family)
            fprintf('[penalty_view] 正式约束族尚未生成，demo 无法选择合格成员。\n');
            close(fig);
            return;
        end
        [i1,i2] = pickTwoMembers(family);

        ps = cell(1,2); sims = cell(1,2);
        for jj = 1:2
            if jj == 1, idx = i1; else, idx = i2; end
            pj = withDefaults(penalty_scenario('defaults'));
            pj.windY = family(idx).windY;
            pj.density = family(idx).density;
            pj.Cd = family(idx).Cd;
            pj.tolerance = p.tolerance;
            ps{jj} = pj;
            sims{jj} = penalty_scenario('start',pj);
            assertSim(sims{jj});
        end

        gifPath = fullfile(outImgDir,'约束轨迹簇演示.gif');
        pngPath = fullfile(outImgDir,'约束轨迹簇对比.png');
        nF = 80; half = nF/2; gifDur = 5.0;
        for jj = 1:2
            selectScenario(sims{jj},ps{jj});
            if jj == 1
                setStatus('约束轨迹簇演示：第 1 组全程参数；环境从出脚起全程生效。');
            else
                setStatus('切换：另一组全程参数（第 2 组）；环境从出脚起全程生效。');
            end
            for kk = 1:half
                tNow = sim.tEnd*(kk-1)/(half-1);
                if kk==29,set(hView,'Value',2);applyView();end
                drawAt(tNow);
                drawnow;
                writeGifFrame(fig,gifPath,(jj-1)*half+kk,gifDur/nF);
            end
        end

        drawComparison(sims);
        setStatus('约束轨迹簇对比：两组差异较大的合格全程参数。');
        drawnow;
        fr = getframe(fig);
        imwrite(frame2im(fr),pngPath);

        G = cell(1,2);
        for jj = 1:2
            G{jj} = buildGroupMeta(sims{jj},ps{jj});
        end
        J = struct();
        J.generatedAt = datestr(now,31);
        J.exportMode = 'demo';
        J.tolerance = p.tolerance;
        J.family = struct('qualified',familyQualified,'drawn',familyDrawn, ...
            'note','等间隔抽样绘制，不是概率，也不是置信区间。');
        J.groups = G;
        J.outputs = struct('gif',gifPath,'png',pngPath);
        J.note = ['两组均为同一观测约束下的合格全程参数；差异来自全程环境/气动参数，' ...
            '不是飞行中途的突然阵风。本文件不覆盖正式导出数据。'];
        writeJson(fullfile(outDataDir,'轨迹簇演示参数.json'),J);
        close(fig);
        fprintf('[penalty_view] demo 完成：%s，%s。\n',gifPath,pngPath);
    end

    function selectScenario(simG,pG)
        sim = simG; p = pG;
        setappdata(fig,'Scenario',sim);
        setappdata(fig,'Parameters',p);
        refreshCurve(); refreshObservations();
        syncControls(); resetTimeZero();
        set(hView,'Value',1);applyStyles(); applyView(); drawAt(0);
    end

    function drawComparison(simsCell)
        for jj = 1:2
            c = penalty_scenario('curve',simsCell{jj});
            set(hCmp(jj),'XData',c.state(:,1),'YData',c.state(:,2),'ZData',c.state(:,3));
            ep = simsCell{jj}.endpoint(:);
            if strcmp(simsCell{jj}.status,'plane')
                set(hCmpGoal(jj),'XData',ep(1),'YData',ep(2),'ZData',ep(3),'Visible','on');
                set(hCmpTxt(jj),'String',sprintf('第%d组交点\n(11, %.3f, %.3f) m',jj,ep(2),ep(3)));
            end
            set(hCmpTxt(jj),'Position',[ep(1)-.7,ep(2)+.08,ep(3)+(.20-.35*(jj-1))]);
        end
        set([hPlayed hFuture hCur hGoalPoint hGoalCross hGoalText],'Visible','off');
        set([hqG hqD hqM],'Visible','off');
        set(hBall,'Visible','off');
        set(hObs,'Visible','on'); set(hFitObs,'Visible','on');
        set(hCmp,'Visible','on'); set(hCmpTxt,'Visible','on');
        set(ax3.Title,'String','约束轨迹簇对比：两组差异较大的合格全程参数');
        set(hView,'Value',2);applyView();view(ax3,-28,20);
    end

    % ================= 小工具 =================
    function b = busyNow()
        b = false;
        try, b = logical(getappdata(fig,'Busy')); catch, end
    end

    function k = paramIndex(name)
        k = find(strcmp(name,PRM.names),1);
        assert(~isempty(k),'未知参数。');
    end
end

% ================= 局部函数 =================
function p = withDefaults(p)
if ~isstruct(p), p = struct(); end
defs = struct('windY',0,'density',1.225,'Cd',0.370687,'tolerance',0.005,'maxStep',0.002);
f = fieldnames(defs);
for k = 1:numel(f)
    if ~isfield(p,f{k}) || ~isnumeric(p.(f{k})) || ~isfinite(p.(f{k}))
        p.(f{k}) = defs.(f{k});
    end
end
end

function assertSim(s)
assert(isstruct(s),'penalty_scenario(''start'') 未返回结构体。');
need = {'params','theta','initialSpeed','fit','status','tEnd','endpoint','base','reference'};
for k = 1:numel(need)
    assert(isfield(s,need{k}),['仿真结果缺少字段：' need{k}]);
end
assert(isfield(s.fit,'feasible') && isfield(s.fit,'maxResidual'), ...
    'sim.fit 缺少 feasible / maxResidual。');
end

function s = fmt(v,dp)
if ~isfinite(v)
    s = 'NaN';
else
    s = sprintf('%.*f',dp,v);
end
end

function txt = statusChinese(st)
switch lower(char(st))
    case 'plane',   txt = '球门平面';
    case 'ground',  txt = '触地';
    case 'timeout', txt = '超时';
    otherwise,      txt = char(st);
end
end

function stopOnly(tm)
if ~isempty(tm) && isvalid(tm)
    try, stop(tm); catch, end
end
end

function stopTimer(tm)
if ~isempty(tm) && isvalid(tm)
    try, stop(tm); catch, end
    try, delete(tm); catch, end
end
end

function writeGifFrame(figH,path,k,delay)
fr = getframe(figH);
[A,map] = rgb2ind(frame2im(fr),256);
if k == 1
    imwrite(A,map,path,'gif','LoopCount',Inf,'DelayTime',delay);
else
    imwrite(A,map,path,'gif','WriteMode','append','DelayTime',delay);
end
end

function writeJson(path,J)
fid = fopen(path,'w','n','UTF-8');
assert(fid >= 0,'无法写入 %s',path);
fprintf(fid,'%s',jsonencode(J,'PrettyPrint',true));
fclose(fid);
end

function J = buildMeta(s,pS,exportMode,outputs)
J = struct();
J.generatedAt = datestr(now,31);
J.exportMode = exportMode;
J.parameters = pS;
J.theta = s.theta(:)';
J.admissibility=s.admissibility;
J.initialSpeed = s.initialSpeed;
J.fit = struct('maxResidual',fitField(s,'maxResidual',NaN), ...
    'feasible',isfield(s.fit,'feasible') && s.fit.feasible, ...
    'residual',fitField(s,'residual',[]), ...
    'observedPositions',fitField(s,'observedPositions',[]), ...
    'exitflag',fitField(s,'exitflag',NaN), ...
    'elapsedSeconds',fitField(s,'elapsedSeconds',NaN), ...
    'kineticEnergy',fitField(s,'kineticEnergy',NaN));
J.endpoint = s.endpoint(:)';
J.status = s.status;
J.tEnd = s.tEnd;
J.isOpening = s.isOpening;
J.referenceEndpoint = s.reference.endpoint(:)';
J.outputs = outputs;
J.note = ['参数、拟合与末端来自 penalty_scenario 的全程联合反演；' ...
    '环境从出脚 t=0 起全程生效。'];
end

function G = buildGroupMeta(s,pS)
G = struct();
G.parameters = pS;
G.theta = s.theta(:)';
G.admissibility=s.admissibility;
G.initialSpeed = s.initialSpeed;
G.kineticEnergy = fitField(s,'kineticEnergy',NaN);
G.fit = struct('maxResidual',fitField(s,'maxResidual',NaN), ...
    'feasible',isfield(s.fit,'feasible') && s.fit.feasible, ...
    'residual',fitField(s,'residual',[]), ...
    'observedPositions',fitField(s,'observedPositions',[]), ...
    'exitflag',fitField(s,'exitflag',NaN), ...
    'elapsedSeconds',fitField(s,'elapsedSeconds',NaN));
G.endpoint = s.endpoint(:)';
G.status = s.status;
G.tEnd = s.tEnd;
G.isOpening = s.isOpening;
end

function v = fitField(s,name,default)
v = default;
if isfield(s.fit,name) && ~isempty(s.fit.(name))
    v = s.fit.(name);
end
end

function fam = loadFamilyDirect(tol)
fam = [];
try
    f = penalty_scenario('family',tol);
    if isstruct(f), fam = f; end
catch
end
if ~isempty(fam) && isfield(fam,'accepted')
    fam = fam([fam.accepted] ~= 0);
end
end

function [i1,i2] = pickTwoMembers(fam)
E = nan(numel(fam),3);
for k = 1:numel(fam)
    e = fam(k).endpoint(:);
    if numel(e) >= 3, E(k,:) = e(1:3)'; end
end
ok = all(isfinite(E),2);
if nnz(ok) < 2
    i1 = 1; i2 = min(2,numel(fam));
    return;
end
Eok = E(ok,:); okIdx = find(ok);
c = mean(Eok,1);
[~,a] = max(sum((Eok-c).^2,2));
i1 = okIdx(a);
d = sum((Eok-Eok(a,:)).^2,2); d(a) = -Inf;
[~,b] = max(d);
i2 = okIdx(b);
end

function runSelftest()
p = penalty_scenario('defaults');
need = {'windY','density','Cd','tolerance','maxStep'};
for k = 1:numel(need)
    assert(isfield(p,need{k}) && isscalar(p.(need{k})),['defaults 缺少参数 ' need{k}]);
end
assert(abs(p.windY) < 1e-12 && abs(p.density-1.225) < 1e-9);
assert(abs(p.tolerance-0.005) < 1e-12 && abs(p.maxStep-0.002) < 1e-12);
s = penalty_scenario('start',p);
assertSim(s);
assert(numel(s.theta) == 6,'theta 必须为 [vx vy vz Cd by bz]。');
c = penalty_scenario('curve',s);
assert(~isempty(c.t) && abs(c.t(1)) < 1e-9,'curve 必须从 t=0 开始。');
assert(abs(c.t(end)-s.tEnd) < 1e-6,'curve 终点必须等于 tEnd。');
sm = penalty_scenario('sample',s,0);
assert(numel(sm.state) == 6 && isfield(sm,'forces') && ...
    all(isfield(sm.forces,{'gravity','drag','magnus'})),'sample 返回字段不完整。');
fam = penalty_scenario('family',p.tolerance);
assert(isempty(fam) || isstruct(fam),'family 必须返回结构体数组或空数组。');
fprintf('[penalty_view] selftest 通过：curve %d 点，theta 六维，family %d 条。\n', ...
    numel(c.t),numel(fam));
end

function text=candidateStatus(s)
if ~s.fit.observationCompatible
 text='该候选未满足观测容差。';
elseif ~s.fit.physicalAdmissible
 text='该候选超出物理筛选范围，不纳入轨迹簇。';
else
 text='满足观测容差与物理筛选。';
end
end
