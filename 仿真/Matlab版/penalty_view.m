function penalty_view(mode)
%PENALTY_VIEW 点球实时阵风/气动参数仿真界面
%   penalty_view                 交互界面（单一可旋转三维图 + 参数面板）
%   penalty_view('export')       导出默认无风 80 帧 GIF / 最终 PNG
%   penalty_view('demo')         导出"相同观测阵风对比" GIF / PNG
%   penalty_view('selftest')     轻量接口自检（不建界面）
%
%   动力学只调用 penalty_scenario 的稳定接口，本文件不实现第二套物理。
%   conditioned 模式前 0.10 s 与基准一致；所有环境/气动改变最早从 0.10 s 生效。
%   mass/radius/mode/speedScale/yawDeg/elevationDeg 只在重开球时生效。

if nargin < 1 || isempty(mode), mode = 'interactive'; end
mode = lower(char(mode));
assert(any(strcmp(mode,{'interactive','export','demo','selftest'})), ...
    'mode 必须是 interactive / export / demo / selftest。');

thisFile = mfilename('fullpath'); outDir = fileparts(thisFile);
outDataDir = fullfile(outDir,'data'); outImgDir = fullfile(outDir,'img');

% ---------- 引擎 ----------
p = penalty_scenario('defaults');
sim = penalty_scenario('start',p);
assert(isstruct(sim) && isfield(sim,'tEnd'), '缺少 penalty_scenario 引擎或返回格式不符。');

if strcmp(mode,'selftest')
    runSelftest(sim);
    return;
end

if any(strcmp(mode,{'export','demo'}))
    if ~exist(outImgDir,'dir'), mkdir(outImgDir); end
    if ~exist(outDataDir,'dir'), mkdir(outDataDir); end
end

% ---------- 基准与观测 ----------
refT = sim.reference.t(:); refS = sim.reference.state;
if size(refS,1) ~= numel(refT), refS = refS.'; end
obsT = sim.base.observations.times(:);
obs  = sim.base.observations.positions;
if size(obs,1) ~= numel(obsT), obs = obs.'; end

% 球门几何（与正式结果一致：x=11 m，半宽 3.66 m，高 2.44 m）
xPlane = 11; goalHW = 3.66; goalH = 2.44;

% ---------- 状态变量 ----------
playing = false; speed = 0.1; t0 = 0; tClock = []; timerH = [];
playbackTime = 0;
trajT = []; trajS = [];
statusText = '就绪'; infoLines = {'','',''};
hS = struct(); hE = struct(); listeners = {};
loMap = struct(); hiMap = struct(); dpMap = struct(); labelMap = struct();
scaleW = 0.05;      % 风箭头：0.05 m / (m/s)
scaleF = 0.25;      % 力箭头：0.25 m / N
fnUI = 'Microsoft YaHei';

% ---------- 图窗与唯一三维轴 ----------
fig = figure('Name','点球实时参数仿真','Color','w','Units','pixels', ...
    'Position',[40 40 1440 860],'NumberTitle','off');
if any(strcmp(mode,{'export','demo'})),set(fig,'Visible','off');end
ax3 = axes(fig,'Units','normalized','Position',[0.405 0.205 0.565 0.750], ...
    'Tag','Trajectory3D');
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on'); view(ax3,-28,16);
set(ax3,'FontName',fnUI);

patch(ax3,[9 11.4 11.4 9],[-0.4 -0.4 3.8 3.8],[0 0 0 0],[0.42 0.72 0.42], ...
    'FaceColor',[0.42 0.72 0.42],'EdgeColor','none','FaceAlpha',0.10, ...
    'HandleVisibility','off');
plot3(ax3,xPlane*[1 1], goalHW*[1 1],[0 goalH],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],-goalHW*[1 1],[0 goalH],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],[-goalHW goalHW],goalH*[1 1],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],[-goalHW goalHW],[0 0],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',1.2,'HandleVisibility','off');

% 基准灰虚线 / 已飞蓝实线 / 预测橙虚线
hRef    = plot3(ax3,refS(:,1),refS(:,2),refS(:,3),'--','Color',[0.55 0.55 0.55], ...
    'LineWidth',1.2,'DisplayName','基准（无风参考）');
hPlayed = plot3(ax3,nan,nan,nan,'-','Color',[0 0.35 0.85],'LineWidth',2.2, ...
    'DisplayName','已飞行');
hFuture = plot3(ax3,nan,nan,nan,'--','Color',[0.95 0.55 0.08],'LineWidth',1.6, ...
    'DisplayName','预测未来');

% 观测点
hObs = plot3(ax3,obs(:,1),obs(:,2),obs(:,3),'ko','MarkerFaceColor','k', ...
    'MarkerSize',6,'DisplayName','观测点');
plot3(ax3,obs(end,1),obs(end,2),obs(end,3),'rs','MarkerSize',9,'LineWidth',1.5, ...
    'HandleVisibility','off');
text(ax3,obs(1,1)+0.10,obs(1,2)-0.85,obs(1,3)+0.55,'观测 t = 0.05 s', ...
    'FontName',fnUI,'FontSize',9,'Clipping','on');
hObsNote = text(ax3,obs(end,1)-0.50,obs(end,2)+1.20,obs(end,3)+0.40, ...
    '观测 t = 0.10 s', 'FontName',fnUI,'FontSize',9,'Clipping','on');

% 足球（真实半径、固定图案、不旋转）
[Xb,Yb,Zb] = sphere(20); ang = atan2(Yb,Xb);
Cb = 0.92*ones(size(Zb)); Cb(mod(floor((ang+pi)/(pi/4)) + (Zb>0),2)==0) = 0.04;
hBall = surf(ax3,Xb,Yb,Zb,Cb,'EdgeColor','none','HandleVisibility','off');
colormap(ax3,gray(64));

% 风箭头（青色，独立单位 m/(m/s)）
hWind = quiver3(ax3,0,0,0,0,0,0,0,'Color',[0 0.72 0.82],'LineWidth',2.0, ...
    'MaxHeadSize',0.8,'HandleVisibility','off');
% 三力箭头（红重力 / 绿阻力 / 蓝 Magnus；0.25 m/N）
hqG = quiver3(ax3,0,0,0,0,0,0,0,'Color',[0.85 0.2 0.2],'LineWidth',1.5,'HandleVisibility','off');
hqD = quiver3(ax3,0,0,0,0,0,0,0,'Color',[0.2 0.65 0.2],'LineWidth',1.5,'HandleVisibility','off');
hqM = quiver3(ax3,0,0,0,0,0,0,0,'Color',[0.2 0.3 0.9],'LineWidth',1.5,'HandleVisibility','off');
set([hqG hqD hqM],'Visible','off');

xlabel(ax3,'x (m)'); ylabel(ax3,'y (m)'); zlabel(ax3,'z (m)');
title(ax3,'点球飞行三维轨迹（灰虚线=基准，蓝实线=已飞，橙虚线=预测）', ...
    'FontName',fnUI,'FontSize',11);
lg = legend(ax3,{'基准（无风参考）','已飞行','预测未来','观测点'}, ...
    'Location','northwest','FontSize',9,'Box','off');
set(lg,'AutoUpdate','off');

% ---------- 左侧参数面板（三个 tab，不新增 axes） ----------
tg = uitabgroup(fig,'Units','normalized','Position',[0.012 0.155 0.325 0.815]);
tabEnv  = uitab(tg,'Title','环境/阵风');
tabAero = uitab(tg,'Title','气动');
tabKick = uitab(tg,'Title','起脚/模式');

makeRow(tabEnv,'windX','顺逆风 X (m/s)',-15,15,1,1);
makeRow(tabEnv,'windY','横风 Y (m/s)',-15,15,1,2);
makeRow(tabEnv,'windZ','升降气流 Z (m/s)',-8,8,2,3);
makeRow(tabEnv,'gustY','阵风 Y 峰值 (m/s)',-20,20,1,4);
makeRow(tabEnv,'gustZ','阵风 Z 峰值 (m/s)',-10,10,2,5);
makeRow(tabEnv,'gustStart','阵风起始 (s)',0.1,0.6,2,6);
makeRow(tabEnv,'gustDuration','阵风持续 (s)',0.02,0.5,2,7);
makeRow(tabEnv,'turbulence','扰动风幅值 (m/s)',0,5,2,8);

makeRow(tabAero,'CdScale','阻力系数倍率',0.3,2,2,1);
makeRow(tabAero,'cdExponent','阻力速度指数',-1,1,2,2);
makeRow(tabAero,'spinScale','有效侧力倍率',-2,2,2,3);
makeRow(tabAero,'bx','有效 Magnus 系数 bx (1/s)',-2,2,2,4);
makeRow(tabAero,'spinTau','侧力衰减时间 (s，0=不衰减)',0,2,2,5);
makeRow(tabAero,'density','空气密度 (kg/m^3)',0.9,1.4,3,6);

% 起脚 tab：模式下拉
y0 = 0.88;
uicontrol(tabKick,'Style','text','Units','normalized','Position',[0.03 y0 0.35 0.06], ...
    'String','飞行模式','HorizontalAlignment','left','FontName',fnUI,'FontSize',9);
hMode = uicontrol(tabKick,'Style','popupmenu','Units','normalized', ...
    'Position',[0.40 y0 0.57 0.06],'Tag','ModeSelect','String',{'固定题给观测','自由起脚实验'}, ...
    'Value',1,'FontName',fnUI,'FontSize',9,'Callback',@(~,~)cbMode());
makeRow(tabKick,'speedScale','初速倍率',0.6,1.4,2,2);
makeRow(tabKick,'yawDeg','偏航角 (°)',-15,15,1,3);
makeRow(tabKick,'elevationDeg','仰角 (°)',-10,10,1,4);
makeRow(tabKick,'mass','球质量 (kg)',0.41,0.45,3,5);
makeRow(tabKick,'radius','球半径 (m)',0.108,0.112,4,6);
uicontrol(tabKick,'Style','text','Units','normalized','Position',[0.03 0.05 0.94 0.22], ...
    'String',['说明：模式、初速倍率、偏航角、仰角、质量、半径只在新开球时生效；' ...
    '改动后自动暂停并重开。固定观测模式下初速/球属性控件禁用，' ...
    '固定观测模式锁定初速和球属性；自由起脚会改变题给观测。'], ...
    'HorizontalAlignment','left','FontName',fnUI,'FontSize',9);

% ---------- 底部控制区 ----------
hTime = uicontrol(fig,'Style','slider','Tag','TimeSlider','Units','normalized', ...
    'Position',[0.020 0.104 0.380 0.032],'Min',0,'Max',max(sim.tEnd,1e-6), ...
    'Value',0,'SliderStep',[0.001 0.05],'Callback',@(~,~)cbScrub());
hPlay = uicontrol(fig,'Style','pushbutton','Tag','PlayButton','Units','normalized', ...
    'Position',[0.410 0.095 0.052 0.045],'String','播放','Callback',@(~,~)cbPlay(), ...
    'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','pushbutton','Tag','ReplayButton','Units','normalized', ...
    'Position',[0.466 0.095 0.052 0.045],'String','重播','Callback',@(~,~)cbReplay(), ...
    'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','popupmenu','Tag','SpeedSelect','Units','normalized', ...
    'Position',[0.522 0.095 0.055 0.045],'String',{'0.05x','0.1x','0.25x','1x'}, ...
    'Value',2,'Callback',@(~,~)cbSpeed(),'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','popupmenu','Tag','ViewSelect','Units','normalized', ...
    'Position',[0.581 0.095 0.115 0.045], ...
    'String',{'全程','末段放大','自动包含未来'},'Value',3, ...
    'Callback',@(~,~)cbView(),'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','checkbox','Tag','ForceToggle','Units','normalized', ...
    'Position',[0.702 0.098 0.070 0.040],'String','受力箭头','Value',0, ...
    'Callback',@(~,~)cbForce(),'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','checkbox','Tag','WindToggle','Units','normalized', ...
    'Position',[0.776 0.098 0.068 0.040],'String','风箭头','Value',1, ...
    'Callback',@(~,~)cbWind(),'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','pushbutton','Tag','ResetButton','Units','normalized', ...
    'Position',[0.848 0.095 0.060 0.045],'String','全部复位', ...
    'Callback',@(~,~)cbReset(),'FontName',fnUI,'FontSize',9);
uicontrol(fig,'Style','pushbutton','Tag','WindDemoButton','Units','normalized', ...
    'Position',[0.912 0.095 0.078 0.045],'String','阵风对比', ...
    'Callback',@(~,~)cbWindDemo(),'FontName',fnUI,'FontSize',9);

hInfo = uicontrol(fig,'Style','text','Units','normalized', ...
    'Position',[0.018 0.004 0.968 0.085],'String','', ...
    'HorizontalAlignment','left','FontName',fnUI,'FontSize',8);

setappdata(fig,'Scenario',sim);
setappdata(fig,'CurrentTime',0);
setappdata(fig,'Parameters',p);
setappdata(fig,'ParamListeners',listeners);
setappdata(fig,'Busy',false);

% ---------- 初始化 ----------
afterRebuild();
set(fig,'CloseRequestFcn',@(~,~)cbClose());

if any(strcmp(mode,{'export','demo'}))
    if strcmp(mode,'demo')
        p = penalty_scenario('defaults');
        p.gustY = 12; p.gustStart = 0.15; p.gustDuration = 0.20;
        sim = penalty_scenario('start',p);
        syncControls(); afterRebuild();
    end
    exportMedia(mode);
    return;
end

% ---------- 交互播放 ----------
timerH = timer('ExecutionMode','fixedSpacing','Period',0.04,'BusyMode','drop', ...
    'TimerFcn',@(~,~)cbTick(),'Name','PenaltyRealtimePlayback');
set(fig,'DeleteFcn',@(~,~)stopTimer(timerH));
rotate3d(fig,'on');
fprintf('[penalty_view] 就绪：T=%.6f s，末端状态 (%s)。\n', sim.tEnd, sim.status);

% ================= 界面构建 =================
    function makeRow(parent,name,label,lo,hi,dp,idx)
        loMap.(name) = lo; hiMap.(name) = hi; dpMap.(name) = dp;
        labelMap.(name) = label;
        y = 0.88 - (idx-1)*0.105;
        uicontrol(parent,'Style','text','Units','normalized','Position',[0.03 y 0.35 0.06], ...
            'String',label,'HorizontalAlignment','left','FontName',fnUI,'FontSize',9);
        hS.(name) = uicontrol(parent,'Style','slider','Units','normalized', ...
            'Position',[0.40 y 0.30 0.055],'Min',lo,'Max',hi,'Value',p.(name), ...
            'Tag',['Param_' name],'SliderStep',[0.01 0.1], ...
            'Callback',@(~,~)cbSlider(name));
        hE.(name) = uicontrol(parent,'Style','edit','Units','normalized', ...
            'Position',[0.73 y 0.24 0.055],'String',fmt(p.(name),dp), ...
            'Tag',['Edit_' name],'Callback',@(~,~)cbEdit(name), ...
            'FontName','Consolas','FontSize',9);
        ls = addlistener(hS.(name),'ContinuousValueChange',@(~,~)cbSliderLive(name));
        listeners{end+1} = ls;
    end

    function afterRebuild()
        p=sim.params;
        setappdata(fig,'Scenario',sim);
        setappdata(fig,'Parameters',p);
        trajT = []; trajS = [];
        refreshCurve();
        set(hTime,'Max',max(sim.tEnd,1e-6),'Value',0);
        playbackTime = 0; setappdata(fig,'CurrentTime',0);
        syncControls();
        if strcmp(p.mode,'free')
            set(hObsNote,'String','观测（自由模式：不再匹配）');
        else
            set(hObsNote,'String','观测 t = 0.10 s（条件模式固定）');
        end
        applyView(get(findobj(fig,'Tag','ViewSelect'),'Value'));
        drawAt(0);
    end

    function refreshCurve()
        c = penalty_scenario('curve',sim);
        t = c.t(:); s = c.state;
        if size(s,1) ~= numel(t), s = s.'; end
        if isempty(t)
            sm = penalty_scenario('sample',sim,0);
            t = 0; s = sm.state(:)';
        elseif t(1) > 1e-9 && ~isempty(trajT)
            k = find(trajT <= t(1)+1e-12);
            t = [trajT(k); t]; s = [trajS(k,:); s];
        end
        [t,iu] = unique(t,'stable'); s = s(iu,:);
        trajT = t; trajS = s;
    end

    function syncControls()
        names = fieldnames(hS);
        for ii = 1:numel(names)
            nm = names{ii};
            set(hS.(nm),'Value',p.(nm));
            set(hE.(nm),'String',fmt(p.(nm),dpMap.(nm)));
        end
        set(hMode,'Value',1+strcmp(p.mode,'free'));
        en = 'on';
        if strcmp(p.mode,'conditioned'), en = 'off'; end
        dirNames = {'speedScale','yawDeg','elevationDeg','mass','radius'};
        for ii = 1:numel(dirNames)
            if isfield(hS,dirNames{ii})
                set([hS.(dirNames{ii}) hE.(dirNames{ii})],'Enable',en);
            end
        end
    end

    function drawAt(tNow)
        tNow = max(0,min(tNow,sim.tEnd));
        playbackTime = tNow; setappdata(fig,'CurrentTime',tNow);
        sm = penalty_scenario('sample',sim,tNow);
        s = sm.state(:)';

        % 已飞/未来曲线
        k = find(trajT <= tNow+1e-12);
        if isempty(k)
            Spast = [s(1:3); s(1:3)];
        else
            Spast = trajS(k,1:3);
            if trajT(k(end)) < tNow-1e-12
                Spast = [Spast(:,1:3); s(1:3)];
            end
        end
        kf = find(trajT > tNow+1e-12);
        Sfut = [s(1:3); trajS(kf,1:3)];
        set(hPlayed,'XData',Spast(:,1),'YData',Spast(:,2),'ZData',Spast(:,3));
        set(hFuture,'XData',Sfut(:,1),'YData',Sfut(:,2),'ZData',Sfut(:,3));

        % 球
        if isfield(sim,'params') && isfield(sim.params,'radius')
            rad = sim.params.radius;
        else
            rad = p.radius;
        end
        set(hBall,'XData',s(1)+rad*Xb,'YData',s(2)+rad*Yb,'ZData',s(3)+rad*Zb);

        % 风箭头（m/s → m，明确比例）
        w = sm.wind(:)';
        set(hWind,'XData',s(1),'YData',s(2),'ZData',s(3), ...
            'UData',scaleW*w(1),'VData',scaleW*w(2),'WData',scaleW*w(3));

        % 三力箭头（N → m）
        F = sm.forces;
        if ~isfield(F,'buoyancy'), F.buoyancy = [0;0;0]; end
        set(hqG,'XData',s(1),'YData',s(2),'ZData',s(3), ...
            'UData',scaleF*F.gravity(1),'VData',scaleF*F.gravity(2),'WData',scaleF*F.gravity(3));
        set(hqD,'XData',s(1),'YData',s(2),'ZData',s(3), ...
            'UData',scaleF*F.drag(1),'VData',scaleF*F.drag(2),'WData',scaleF*F.drag(3));
        set(hqM,'XData',s(1),'YData',s(2),'ZData',s(3), ...
            'UData',scaleF*F.magnus(1),'VData',scaleF*F.magnus(2),'WData',scaleF*F.magnus(3));

        set(hTime,'Value',tNow);

        % 预计终点说明
        ep = sim.endpoint(:)';
        switch sim.status
            case 'plane',  eTxt = sprintf('球门平面 (y,z) = (%.2f, %.2f) m',ep(2),ep(3));
            case 'ground', eTxt = sprintf('触地于 (x,y) = (%.2f, %.2f) m',ep(1),ep(2));
            otherwise,     eTxt = '未到达（时限内未越门/未触地）';
        end
        refEp = sim.reference.endpoint(:)';
        if numel(refEp) >= 3 && all(isfinite(refEp(1:3))) && numel(ep) >= 3 && all(isfinite(ep(1:3)))
            dOff = norm(ep(1:3)-refEp(1:3));
            offTxt = sprintf('%.3f m',dOff);
        else
            offTxt = '—';
        end
        if strcmp(p.mode,'free')
            modeTxt = '自由模式（观测不再匹配）';
        else
            modeTxt = '条件模式（前 0.10 s 沿基准）';
        end
        infoLines = { ...
            sprintf('t = %.3f s ｜ 位置 (%.2f, %.2f, %.2f) m ｜ 速率 %.2f m/s', ...
                tNow,s(1),s(2),s(3),norm(s(4:6))), ...
            sprintf('风 (%+.2f, %+.2f, %+.2f) m/s（|w| = %.2f） ｜ 有效 Cd = %.3f ｜ %s', ...
                w(1),w(2),w(3),norm(w),sm.Cd,modeTxt), ...
            sprintf('三力大小 (N)：重力 %.2f / 阻力 %.2f / 侧力 %.2f',norm(F.gravity),norm(F.drag),norm(F.magnus)), ...
            sprintf('预计：%s ｜ 与基准终点偏移 %s ｜ 风箭头 %.2f m/(m/s)，力箭头 %.2f m/N', ...
                eTxt,offTxt,scaleW,scaleF)};
        renderInfo();
    end

    function renderInfo()
        set(hInfo,'String',[infoLines, {['状态：' statusText]}]);
    end

    function status(msg)
        statusText = msg; renderInfo();
    end

    function applyView(idx)
        if idx == 1
            xr = [-0.6 11.8]; yr = [-4.2 4.2]; zr = [0 3.2];
        elseif idx == 2
            xr = [7 11.6]; yr = [1.2 3.7]; zr = [1.1 2.9];
        else
            X = [refS(:,1); trajS(:,1); 0; 11];
            Y = [refS(:,2); trajS(:,2); -0.3; 0.3];
            Z = [refS(:,3); trajS(:,3); 0];
            xr = [min(X)-0.5, max(X)+0.5];
            yr = [min(Y)-0.4, max(Y)+0.4];
            zr = [min(Z)-0.2, max(Z)+0.4];
            if diff(yr) < 1, yr = yr + [-0.5 0.5]; end
            if diff(zr) < 1, zr = zr + [-0.5 0.5]; end
        end
        xlim(ax3,xr); ylim(ax3,yr); zlim(ax3,zr);
        pbaspect(ax3,[diff(xr) diff(yr) diff(zr)]);   % 单位等比
    end

    % ---------- 参数回调 ----------
    function cbSlider(name)
        if getBusy(fig), return; end
        applyParam(name,get(hS.(name),'Value'));
    end

    function cbSliderLive(name)
        if getBusy(fig), return; end
        applyParam(name,get(hS.(name),'Value'));
    end

    function cbEdit(name)
        if getBusy(fig), return; end
        v = str2double(get(hE.(name),'String'));
        if isnan(v) || ~isfinite(v)
            set(hE.(name),'String',fmt(p.(name),dpMap.(name)));
            status(sprintf('%s 输入无效（非数字/NaN），已恢复旧值。',labelMap.(name)));
            return;
        end
        applyParam(name,v);
    end

    function applyParam(name,val)
        if ~isfield(p,name) || ~isnumeric(val) || ~isscalar(val) || ~isfinite(val)
            status('输入无效，已忽略。'); return;
        end
        if getBusy(fig), return; end
        setappdata(fig,'Busy',true);
        cleanupObj = onCleanup(@() setappdata(fig,'Busy',false)); %#ok<NASGU>

        val = round(val*10^dpMap.(name))/10^dpMap.(name);
        val = min(max(val,loMap.(name)),hiMap.(name));
        set(hS.(name),'Value',val);
        set(hE.(name),'String',fmt(val,dpMap.(name)));

        old = p.(name); p.(name) = val;
        setappdata(fig,'Parameters',p);

        if isRestartField(name)
            stopPlayback();
            sim = penalty_scenario('start',p);
            afterRebuild();
            status(sprintf('%s 已更改：自动暂停并重开球（仅新开球生效）。',labelMap.(name)));
            return;
        end

        if playbackTime >= sim.tEnd - 1e-9
            status(sprintf('%s 已保存；轨迹位于终点，未逆写终态。重播或拖回过去后再改参生效。', ...
                labelMap.(name)));
            return;
        end

        try
            sim = penalty_scenario('branch',sim,playbackTime,p);
        catch err
            p.(name) = old; setappdata(fig,'Parameters',p); syncControls();
            status(['分支求解失败，已恢复旧值：' err.message]);
            return;
        end
        setappdata(fig,'Scenario',sim);
        refreshCurve();
        set(hTime,'Max',max(sim.tEnd,1e-6));
        if playing, t0 = playbackTime; tClock = tic; end   % 重置时钟基点
        drawAt(playbackTime);
        effectiveTime=playbackTime;if strcmp(p.mode,'conditioned'),effectiveTime=max(.1,effectiveTime);end
        status(sprintf('%s = %s 已生效，未来轨迹从 t=%.3f s 重算。', ...
            labelMap.(name),fmt(val,dpMap.(name)),effectiveTime));
    end

    function cbMode()
        modes={'conditioned','free'}; nm=modes{get(hMode,'Value')};
        if strcmp(nm,p.mode) || getBusy(fig), return; end
        setappdata(fig,'Busy',true);
        cleanupObj = onCleanup(@() setappdata(fig,'Busy',false)); %#ok<NASGU>
        p.mode = nm;
        stopPlayback();
        sim = penalty_scenario('start',p);
        afterRebuild();
        status(sprintf('模式已切换为 %s：自动暂停并重开，观测匹配说明已更新。',nm));
    end

    % ---------- 播放 ----------
    function cbPlay()
        if playing
            stopPlayback(); return;
        end
        if playbackTime >= sim.tEnd - 1e-9
            sim=penalty_scenario('start',p); afterRebuild();
        end
        t0 = playbackTime; tClock = tic; playing = true;
        start(timerH); set(hPlay,'String','暂停');
        status('播放中（滑块可连续拖动；改参只重算未来）。');
    end

    function cbTick()
        if ~playing, return; end
        if isempty(tClock), tClock = tic; end
        tn = t0 + speed*toc(tClock);
        if tn >= sim.tEnd
            drawAt(sim.tEnd);
            playing = false; stopOnly(timerH); set(hPlay,'String','播放');
            status('播放结束。可重播，或拖回过去改参数产生新分支。');
            return;
        end
        drawAt(tn);
    end

    function cbReplay()
        stopPlayback();
        sim = penalty_scenario('start',p);   % 保留当前参数，清除手动历史
        afterRebuild();
        status('已重播：保留当前参数重开球，手动历史分支已清除。');
    end

    function cbScrub()
        if playing, stopPlayback(); end
        tn = get(hTime,'Value');
        drawAt(tn);
        status(sprintf('已定位到 t = %.3f s（拖动不改轨迹，可在此改参产生新分支）。',tn));
    end

    function cbSpeed()
        vv = [0.05 0.1 0.25 1];
        speed = vv(get(findobj(fig,'Tag','SpeedSelect'),'Value'));
        t0 = playbackTime; tClock = tic;
    end

    function cbView()
        applyView(get(findobj(fig,'Tag','ViewSelect'),'Value'));
    end

    function cbForce()
        onoff = 'off';
        if get(findobj(fig,'Tag','ForceToggle'),'Value'), onoff = 'on'; end
        set([hqG hqD hqM],'Visible',onoff);
    end

    function cbWind()
        onoff = 'off';
        if get(findobj(fig,'Tag','WindToggle'),'Value'), onoff = 'on'; end
        set(hWind,'Visible',onoff);
    end

    function cbReset()
        stopPlayback();
        p = penalty_scenario('defaults');
        sim = penalty_scenario('start',p);
        afterRebuild();
        status('已全部复位到默认参数并重开球。');
    end

    function cbWindDemo()
        stopPlayback();
        p = penalty_scenario('defaults');
        p.gustY = 12; p.gustStart = 0.15; p.gustDuration = 0.20;
        sim = penalty_scenario('start',p);
        afterRebuild();
        status('阵风对比：gustY=12 m/s，起始 0.15 s，持续 0.20 s；前 0.10 s 与基准相同。');
        cbPlay();
    end

    function stopPlayback()
        playing = false;
        stopOnly(timerH);
        if ishandle(hPlay), set(hPlay,'String','播放'); end
    end

    % ---------- 导出 ----------
    function exportMedia(exportMode)
        if strcmp(exportMode,'demo')
            gifPath = fullfile(outImgDir,'阵风实时仿真.gif');
            pngPath = fullfile(outImgDir,'阵风轨迹对比.png');
        else
            gifPath = fullfile(outImgDir,'点球物理仿真.gif');
            pngPath = fullfile(outImgDir,'点球仿真全景.png');
        end
        nF = 80; gifDur = 4; tMax = sim.tEnd;
        tmpDir = [tempname '_Frames']; mkdir(tmpDir);
        editTime=tMax*35/(nF-1);
        for k = 1:nF
            if strcmp(exportMode,'demo') && k==36
                % 演示真实的中途参数修改：仅改变一个横风分量，不回写历史。
                p.windY=15;sim=penalty_scenario('branch',sim,editTime,p);
                setappdata(fig,'Scenario',sim);setappdata(fig,'Parameters',p);
                refreshCurve();syncControls();set(hTime,'Max',sim.tEnd);
                applyView(get(findobj(fig,'Tag','ViewSelect'),'Value'));
                status('演示：在当前时刻把横风Y改为15 m/s，已飞行状态保持不变。');
            end
            if strcmp(exportMode,'demo') && k>=36
                frameTime=editTime+(sim.tEnd-editTime)*(k-36)/(nF-36);
            else
                frameTime=tMax*(k-1)/(nF-1);
            end
            drawAt(frameTime); drawnow;
            [A,map] = captureRGB(fig,tmpDir,k);
            if k == 1
                imwrite(A,map,gifPath,'gif','LoopCount',Inf,'DelayTime',gifDur/nF);
            else
                imwrite(A,map,gifPath,'gif','WriteMode','append','DelayTime',gifDur/nF);
            end
        end
        drawAt(sim.tEnd); drawnow;
        fr = getframe(fig); imwrite(frame2im(fr),pngPath);
        runName='实时参数导出.json';if strcmp(exportMode,'demo'),runName='阵风参数导出.json';end
        writeRunJson(fullfile(outDataDir,runName),sim,p,exportMode, ...
            struct('gif',gifPath,'png',pngPath));
        if exist(tmpDir,'dir'), rmdir(tmpDir,'s'); end
        close(fig);
        fprintf('[penalty_view] %s 导出完成：T=%.6f s，%s。\n',exportMode,sim.tEnd,sim.status);
    end

    function cbClose()
        try, stopTimer(timerH); catch, end
        try, delete(fig); catch, end
    end
end

% ================= 局部函数 =================
function tf = isRestartField(name)
tf = any(strcmp(name,{'mode','speedScale','yawDeg','elevationDeg','mass','radius'}));
end

function b = getBusy(figH)
b = false;
try
    b = logical(getappdata(figH,'Busy'));
catch
end
end

function s = fmt(v,dp)
if ~isfinite(v), s = 'NaN'; else, s = sprintf('%.*f',dp,v); end
end

function [A,map] = captureRGB(fig,tmpDir,k)
try
    [A,map] = rgb2ind(frame2im(getframe(fig)),256);
catch
    fn = fullfile(tmpDir,sprintf('PenaltyView_tmp_%03d.png',k));
    exportgraphics(fig,fn,'Resolution',120);
    [A,map] = rgb2ind(imread(fn),256); delete(fn);
end
end

function stopOnly(tm)
if ~isempty(tm) && isvalid(tm), try, stop(tm); catch, end, end
end

function stopTimer(tm)
if ~isempty(tm) && isvalid(tm)
    try, stop(tm); catch, end
    try, delete(tm); catch, end
end
end

function writeRunJson(path,sim,p,exportMode,outputs)
J = struct();
J.generatedAt = datestr(now,31);
J.exportMode = exportMode;
J.params = p;
J.status = sim.status;
J.endpoint = struct('t_s',sim.tEnd,'state',sim.endpoint(:)');
J.edits = sim.edits;
J.referenceEndpoint = struct('t_s',sim.reference.tEnd,'state',sim.reference.endpoint(:)');
J.observations = struct('times_s',sim.base.observations.times, ...
    'positions_m',sim.base.observations.positions);
J.outputs = outputs;
J.units = '位置 m；速度 m/s；时间 s；力 N；风箭头 0.05 m/(m/s)；力箭头 0.25 m/N';
J.note = ['参数与事件来自 penalty_scenario；历史力按当时参数采样。' ...
    'b 为有效 Magnus 系数，不代表已测角速度；球皮图案不表示真实自旋。'];
fid = fopen(path,'w','n','UTF-8');
assert(fid >= 0,'无法写入 %s',path);
fprintf(fid,'%s',jsonencode(J,'PrettyPrint',true)); fclose(fid);
end

function runSelftest(sim)
c = penalty_scenario('curve',sim);
assert(~isempty(c.t) && abs(c.t(1)) < 1e-9, 'curve 必须从 0 开始。');
assert(abs(c.t(end)-sim.tEnd) < 1e-6, 'curve 终点与 tEnd 不一致。');
sm = penalty_scenario('sample',sim,0);
assert(numel(sm.state) == 6 && isfield(sm,'forces') && all(isfield(sm.forces, ...
    {'gravity','drag','magnus','buoyancy'})), 'sample 返回字段不完整。');
fprintf('[penalty_view] selftest 通过：curve %d 点，状态六维，力字段齐全。\n',numel(c.t));
end
