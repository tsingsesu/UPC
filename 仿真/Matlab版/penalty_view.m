function penalty_view(mode)
%PENALTY_VIEW 点球动态仿真：交互播放 / 导出（不重拟合，不重跑 Bootstrap）
%   penalty_view / penalty_view('export')，参数读自 罚球结果.json，
%   物理调用 ball_rhs；b=[bx;by;bz] 为有效 Magnus 系数（1/s），不谈真实自旋。
if nargin < 1 || isempty(mode), mode = 'interactive'; end
mode = lower(char(mode));
assert(any(strcmp(mode,{'interactive','export'})), 'mode 必须是 ''interactive'' 或 ''export''。');
isExport = strcmp(mode,'export');
thisFile = mfilename('fullpath'); rootDir = fileparts(fileparts(fileparts(thisFile)));
% 读取：正式计算结果仍在项目根 data/派生数据；导出：演示物写本目录 data/ 与 img/
dataDir = fullfile(rootDir,'data','派生数据'); imgDir = fullfile(rootDir,'img');
outDataDir = fullfile(fileparts(thisFile),'data'); outImgDir = fullfile(fileparts(thisFile),'img');
if isExport
    if ~exist(outImgDir,'dir'), mkdir(outImgDir); end
    if ~exist(outDataDir,'dir'), mkdir(outDataDir); end
end
srcJson = fullfile(dataDir,'罚球结果.json');
assert(exist(srcJson,'file')==2, '缺少已有结果文件：%s', srcJson);
if isExport, assert(usejava('jvm'), 'export 模式需要图形环境（JVM）。'); end
res = jsondecode(fileread(srcJson));

% ---------- 配置（取自已有 JSON） ----------
cfg = struct('R',res.constants.R,'m',res.constants.m,'rho',res.constants.rho, ...
    'g',res.constants.g(:),'r0',res.constants.r0(:),'K',res.constants.K, ...
    'w',[0;0;0],'bx',0,'cdExponent',0,'cdRefSpeed',30,'bDecayTau',inf);
ns = res.numericSettings;
cfg.relTol = ns.predictionRelTol; cfg.absTol = ns.predictionAbsTol;
cfg.odeMaxStep = ns.predictionMaxStep_s; Tmax = ns.maxFlightTime_s;
R = cfg.R; xPlane = ns.goalPlane_x_m;
goalHW = res.constants.goalHalfWidth; goalH = res.constants.goalHeight;
theta = res.fit.theta(:)'; obsT = res.observations.times(:)'; obs = res.observations.positions;
bn = {res.baselines.name};
thetaG = [res.baselines(find(contains(bn,'重力'),1)).v0(:)', 0, 0, 0];
thetaD = res.baselines(find(contains(bn,'仅阻力'),1)).theta(:)';

% ---------- 名义预测、事件与验证（ode45 0→2 s，deval 取 t_eval） ----------
pNom = predictPlane(theta,cfg,xPlane,Tmax);
assert(strcmp(pNom.status,'plane'), '名义轨迹未到达球门平面。');
tAll = unique([linspace(0,pNom.tCross,220), obsT(:)', pNom.tCross])';
SAll = deval(pNom.sol,tAll)'; Pobs = deval(pNom.sol,obsT)';
dObs = max(sqrt(sum((Pobs(:,1:3)-obs).^2,2)));
ePos = norm([pNom.y;pNom.z] - [res.prediction.yGoal;res.prediction.zGoal]);
eT = abs(pNom.tCross - res.prediction.tCross);
assert(dObs < 1e-7, '观测位置误差 %.3g m 超过 1e-7 m。', dObs);
assert(ePos < 1e-7 && eT < 1e-7, '终点与已有 JSON 不一致（%.3g m / %.3g s）。', ePos, eT);
pG = predictPlane(thetaG,cfg,xPlane,Tmax); pD = predictPlane(thetaD,cfg,xPlane,Tmax);
SG = deval(pG.sol,linspace(0,tEndOf(pG,Tmax),160))'; SD = deval(pD.sol,linspace(0,tEndOf(pD,Tmax),160))';

% ---------- 界面（三维 + 俯视双视图，export 时不可见） ----------
fig = figure('Name','点球动态仿真','Color','w','Units','pixels','Position',[40 40 1440 860]);
if isExport, set(fig,'Visible','off'); end
ax3 = axes(fig,'Position',[0.045 0.42 0.44 0.52],'Tag','Trajectory3D');
axTop = axes(fig,'Position',[0.515 0.42 0.44 0.52],'Tag','TrajectoryTop');
hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on'); view(ax3,-28,16);
hold(axTop,'on'); grid(axTop,'on'); box(axTop,'on');
set([ax3 axTop],'FontName','Microsoft YaHei');

% 草坪：仅球门附近小范围淡色
patch(ax3,[9 11.4 11.4 9],[-0.4 -0.4 3.8 3.8],[0 0 0 0],[0.42 0.72 0.42], ...
    'FaceColor',[0.42 0.72 0.42],'EdgeColor','none','FaceAlpha',0.10, ...
    'HandleVisibility','off');
patch(axTop,[9 11.4 11.4 9],[-0.4 -0.4 3.8 3.8],[0.42 0.72 0.42], ...
    'EdgeColor','none','FaceAlpha',0.10,'HandleVisibility','off');

% 球门：只画落在视野内的相关部分（右立柱、横梁、门线）
plot3(ax3,xPlane*[1 1],goalHW*[1 1],[0 goalH],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],[-0.4 goalHW],goalH*[1 1],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',2,'HandleVisibility','off');
plot3(ax3,xPlane*[1 1],[-0.4 goalHW],[0 0],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',1.2,'HandleVisibility','off');
plot(axTop,xPlane*[1 1],[-0.4 goalHW],'-','Color',[0.15 0.15 0.15], ...
    'LineWidth',1.5,'HandleVisibility','off');

% 完整 Magnus 轨迹：已飞行实线 / 未来浅蓝虚线（hFuture 初始即完整轨迹）
hFuture = plot3(ax3,SAll(:,1),SAll(:,2),SAll(:,3),'--','Color',[0.55 0.72 0.92], ...
    'LineWidth',1.2,'DisplayName','未来轨迹');
hPlayed = plot3(ax3,SAll(1,1),SAll(1,2),SAll(1,3),'-','Color',[0 0.35 0.85], ...
    'LineWidth',2.2,'DisplayName','已飞行轨迹');
hG = plot3(ax3,SG(:,1),SG(:,2),SG(:,3),'-.','Color',[0.55 0.55 0.55], ...
    'LineWidth',1,'DisplayName','仅重力基线','Visible','off');
hD = plot3(ax3,SD(:,1),SD(:,2),SD(:,3),'--','Color',[0.85 0.55 0.1], ...
    'LineWidth',1,'DisplayName','仅阻力基线','Visible','off');
hFutureTop = plot(axTop,SAll(:,1),SAll(:,2),'--','Color',[0.55 0.72 0.92], ...
    'LineWidth',1.2,'HandleVisibility','off');
hPlayedTop = plot(axTop,SAll(1,1),SAll(1,2),'-','Color',[0 0.35 0.85], ...
    'LineWidth',2.2,'HandleVisibility','off');
hGTop = plot(axTop,SG(:,1),SG(:,2),'-.','Color',[0.55 0.55 0.55], ...
    'LineWidth',1,'Visible','off','HandleVisibility','off');
hDTop = plot(axTop,SD(:,1),SD(:,2),'--','Color',[0.85 0.55 0.1], ...
    'LineWidth',1,'Visible','off','HandleVisibility','off');

% 观测点（两处文字偏移不同，避免互相遮挡）
hObs = plot3(ax3,obs(:,1),obs(:,2),obs(:,3),'ko','MarkerFaceColor','k', ...
    'MarkerSize',6,'DisplayName','观测点');
plot3(ax3,obs(2,1),obs(2,2),obs(2,3),'rs','MarkerSize',9,'LineWidth',1.5, ...
    'HandleVisibility','off');
plot(axTop,obs(:,1),obs(:,2),'ko','MarkerFaceColor','k','MarkerSize',6, ...
    'HandleVisibility','off');
text(ax3,obs(1,1)+0.10,obs(1,2)-0.85,obs(1,3)+0.55,'观测 t=0.05 s', ...
    'FontName','Microsoft YaHei','FontSize',9,'Clipping','on');
text(ax3,obs(2,1)-0.45,obs(2,2)+1.20,obs(2,3)+0.40,'观测 t=0.10 s（预测起点）', ...
    'FontName','Microsoft YaHei','FontSize',9,'Clipping','on');
% 足球：真实半径 R、固定黑白图案、不旋转（不虚构角速度）
[Xb,Yb,Zb] = sphere(20); ang = atan2(Yb,Xb);
Cb = 0.92*ones(size(Zb)); Cb(mod(floor((ang+pi)/(pi/4)) + (Zb>0),2)==0) = 0.04;
hBall = surf(ax3,cfg.r0(1)+R*Xb,cfg.r0(2)+R*Yb,cfg.r0(3)+R*Zb,Cb,'EdgeColor','none', ...
    'HandleVisibility','off'); colormap(ax3,gray(64));
hCenterLegend = plot3(ax3,nan,nan,nan,'o','MarkerSize',6,'MarkerEdgeColor','k','MarkerFaceColor','w');
hBallTop = plot(axTop,cfg.r0(1),cfg.r0(2),'o','MarkerSize',6, ...
    'MarkerFaceColor',[0 0.35 0.85],'MarkerEdgeColor','k','HandleVisibility','off');

% 受力箭头（默认关闭；开启后按 0.25 m/N 固定比例示意）
sc = 0.25;
hqG = quiver3(ax3,0,0,0,0,0,0,0,'Color',[0.85 0.2 0.2],'LineWidth',1.6,'HandleVisibility','off');
hqD = quiver3(ax3,0,0,0,0,0,0,0,'Color',[0.2 0.65 0.2],'LineWidth',1.6,'HandleVisibility','off');
hqM = quiver3(ax3,0,0,0,0,0,0,0,'Color',[0.2 0.3 0.9],'LineWidth',1.6,'HandleVisibility','off');
set([hqG hqD hqM],'Visible','off');

% 俯视图：门线 x=11 与名义到达点
plot(axTop,xPlane*[1 1],[-0.4 3.8],':','Color',[0.2 0.2 0.2],'LineWidth',1.2, ...
    'HandleVisibility','off');
plot(axTop,xPlane,res.prediction.yGoal,'p','MarkerSize',12,'MarkerFaceColor','y', ...
    'MarkerEdgeColor','k','HandleVisibility','off');
text(axTop,xPlane-3.6,res.prediction.yGoal+0.28,sprintf('门线 x=11 m；名义到达 y=%.2f m', ...
    res.prediction.yGoal),'FontName','Microsoft YaHei','FontSize',9);

xlabel(ax3,'x (m)'); ylabel(ax3,'y (m)'); zlabel(ax3,'z (m)');
title(ax3,'点球飞行三维轨迹（实线已飞行，虚线模型推算）','FontName','Microsoft YaHei');
xlabel(axTop,'x (m)'); ylabel(axTop,'y (m)');
title(axTop,'俯视图（x-y，等比例）','FontName','Microsoft YaHei');

% 共享图例：放在两图下方、控件上方；只列 6 项
lg = legend(ax3,[hPlayed hFuture hG hD hCenterLegend hObs], ...
    {'已飞行轨迹','未来轨迹（推算）','仅重力基线','仅阻力基线','球心','观测点'}, ...
    'Orientation','horizontal','NumColumns',6,'FontSize',9);
set(lg,'Units','normalized','Position',[0.045 0.335 0.91 0.05], ...
    'AutoUpdate','off','Box','off');

% ---------- 控件（底部控制区，不遮挡 axes） ----------
hSlider = uicontrol(fig,'Style','slider','Tag','TimeSlider','Units','normalized', ...
    'Position',[0.035 0.245 0.345 0.035],'Min',0,'Max',pNom.tCross,'Value',0, ...
    'Callback',@cbSlider);
hPlay = uicontrol(fig,'Style','pushbutton','Tag','PlayButton','Units','normalized', ...
    'Position',[0.390 0.238 0.060 0.045],'String','播放','Callback',@cbPlay, ...
    'FontName','Microsoft YaHei','FontSize',9);
uicontrol(fig,'Style','pushbutton','Tag','ReplayButton','Units','normalized', ...
    'Position',[0.456 0.238 0.060 0.045],'String','重播','Callback',@cbReplay, ...
    'FontName','Microsoft YaHei','FontSize',9);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.528 0.283 0.075 0.028], ...
    'String','播放速度','FontName','Microsoft YaHei','FontSize',9);
uicontrol(fig,'Style','popupmenu','Tag','SpeedSelect','Units','normalized', ...
    'Position',[0.528 0.238 0.075 0.045],'String',{'0.05x','0.1x','0.25x','1x'}, ...
    'Value',2,'Callback',@cbSpeed,'FontName','Microsoft YaHei','FontSize',9);
uicontrol(fig,'Style','text','Units','normalized','Position',[0.615 0.283 0.070 0.028], ...
    'String','视野','FontName','Microsoft YaHei','FontSize',9);
uicontrol(fig,'Style','popupmenu','Tag','ViewSelect','Units','normalized', ...
    'Position',[0.615 0.238 0.135 0.045], ...
    'String',{'全程（x -0.4~11.6 m）','末段放大（x 7~11.6 m）'},'Value',1, ...
    'Callback',@cbView,'FontName','Microsoft YaHei','FontSize',9);
hView = findobj(fig,'Tag','ViewSelect');
hCmp = uicontrol(fig,'Style','checkbox','Tag','BaselineToggle','Units','normalized', ...
    'Position',[0.765 0.240 0.095 0.040],'String','叠加基线','Value',0, ...
    'Callback',@cbCmp,'FontName','Microsoft YaHei','FontSize',9);
hForce = uicontrol(fig,'Style','checkbox','Tag','ForceToggle','Units','normalized', ...
    'Position',[0.875 0.240 0.110 0.040],'String','受力箭头','Value',0, ...
    'Callback',@cbForce,'FontName','Microsoft YaHei','FontSize',9);
hTL = uicontrol(fig,'Style','text','Units','normalized','Position',[0.035 0.175 0.94 0.035], ...
    'String','t = 0 s','FontName','Microsoft YaHei','FontSize',9,'HorizontalAlignment','left');
hFL = uicontrol(fig,'Style','text','Units','normalized','Position',[0.035 0.128 0.94 0.035], ...
    'String','力 (N)','FontName','Microsoft YaHei','FontSize',9,'HorizontalAlignment','left');
uicontrol(fig,'Style','text','Units','normalized','Position',[0.035 0.015 0.94 0.085], ...
    'String',['早期观测之后为模型推算；b 是有效 Magnus 系数，球皮图案不表示真实自旋。' ...
    ' 视野下拉只改变显示范围，不改变数据；受力箭头默认关闭，开启后按 0.25 m/N 固定比例仅作示意。'], ...
    'FontName','Microsoft YaHei','FontSize',9,'HorizontalAlignment','left');
applyView(1);
playing = false; speed = 0.1; t0 = 0; tClock = tic; timerH = [];
set(fig,'CloseRequestFcn',@cbClose); updateFrame(0);
if ~isExport
    timerH = timer('ExecutionMode','fixedSpacing','Period',0.04,'BusyMode','drop', ...
        'TimerFcn',@cbTick,'Name','PenaltyViewPlayback');
    set(fig,'DeleteFcn',@(~,~)stopTimer(timerH));
    fprintf('[penalty_view] 交互界面就绪：T=%.9f s，y=%.6f m，z=%.6f m。\n', ...
        pNom.tCross,pNom.y,pNom.z);
    return;
end

% ---------- export：GIF、最终 PNG、时序 CSV、验证 JSON ----------
gifPath = fullfile(outImgDir,'点球物理仿真.gif'); pngPath = fullfile(outImgDir,'点球仿真全景.png');
csvPath = fullfile(outDataDir,'点球仿真时序.csv'); verPath = fullfile(outDataDir,'点球仿真验证.json');
set(hCmp,'Value',1); cbCmp(hCmp);               % 导出：基线开启
set(hForce,'Value',0); cbForce(hForce);         % 导出：箭头关闭
set(hView,'Value',1); cbView(hView);             % 导出：全程
nF = 80; gifDur = 4; tmpDir = [tempname '_Frames']; mkdir(tmpDir);
for k = 1:nF
    updateFrame(tAll(end)*(k-1)/(nF-1)); drawnow;
    [A,map] = captureRGB(fig,tmpDir,k);
    if k == 1
        imwrite(A,map,gifPath,'gif','LoopCount',Inf,'DelayTime',gifDur/nF);
    else
        imwrite(A,map,gifPath,'gif','WriteMode','append','DelayTime',gifDur/nF);
    end
end
updateFrame(tAll(end)); drawnow; frFinal=getframe(fig); imwrite(frame2im(frFinal),pngPath);
tf = linspace(0,tAll(end),nF)'; Sf = deval(pNom.sol,tf)';
[Fg,Fd,Fm] = forceSeries(tf,Sf,theta,cfg);
Tcsv = table(tf,Sf(:,1),Sf(:,2),Sf(:,3),Sf(:,4),Sf(:,5),Sf(:,6), ...
    Fg(:,1),Fg(:,2),Fg(:,3),Fd(:,1),Fd(:,2),Fd(:,3),Fm(:,1),Fm(:,2),Fm(:,3), ...
    'VariableNames',{'t_s','x_m','y_m','z_m','vx_mps','vy_mps','vz_mps','Fgx_N', ...
    'Fgy_N','Fgz_N','Fdx_N','Fdy_N','Fdz_N','Fmx_N','Fmy_N','Fmz_N'});
writetable(Tcsv,csvPath);
ver = struct();
ver.status = pNom.status; ver.integrator = ns.predictionIntegrator;
ver.relTol = cfg.relTol; ver.absTol = cfg.absTol; ver.maxStep_s = cfg.odeMaxStep;
ver.eventRule = ['ode45 从 0 至 2 s；事件 x=11 m（方向+1，终止）；z=R（方向-1，终止）；' ...
    't<1e-7 s 时地面事件值取为正数，避免起步初始接触误终止'];
ver.events = struct('times_s',pNom.te,'codes',pNom.ie,'codesMeaning','1=x=11 m 球门平面；2=触地 z=R');
ver.observationMaxError_m = dObs;
ver.observed_m = obs; ver.predictedAtObs_m = Pobs(:,1:3);
ver.endpointPositionError_m = ePos; ver.endpointTimeError_s = eT;
ver.endpoint_m = struct('x',pNom.ye(1,1),'y',pNom.y,'z',pNom.z,'t_s',pNom.tCross);
ver.existingJSON = struct('yGoal',res.prediction.yGoal,'zGoal',res.prediction.zGoal, ...
    'tCross_s',res.prediction.tCross);
ver.theta = theta;
ver.effectiveMagnusCoefficient = struct('bx',cfg.bx,'by',theta(5),'bz',theta(6),'unit','1/s');
ver.baselines = struct('gravityV0',thetaG(1:3),'dragV0Cd',thetaD([1 2 3 4]));
ver.forceArrowScale_mPerN = sc;
ver.viewLayout = '3D+XY'; ver.predictionPanel = false;
ver.viewRanges = struct('full',struct('x',[-0.4 11.6],'y',[-0.4 3.8],'z',[0 3]), ...
    'zoom',struct('x',[7 11.6],'y',[1.3 3.7],'z',[1.2 2.8]));
ver.viewUsedExport = 'full';
ver.units = '位置 m；速度 m/s；时间 s；力 N';
ver.frames = nF; ver.gifNominalDuration_s = gifDur; ver.gifDelay_s = gifDur/nF;
ver.outputs = struct('gif',gifPath,'png',pngPath,'csv',csvPath,'json',verPath);
ver.note = ['只显示名义有效 Magnus 模型；观测之外为模型推算；球皮图案不表示真实自旋，' ...
    '未虚构角速度；力矢量仅按 0.25 m/N 示意。'];
fid = fopen(verPath,'w','n','UTF-8'); assert(fid >= 0,'无法写入验证 JSON。');
fprintf(fid,'%s',jsonencode(ver,'PrettyPrint',true)); fclose(fid);
if exist(tmpDir,'dir'), rmdir(tmpDir,'s'); end
close(fig);
fprintf('[penalty_view] export 完成：观测误差 %.3g m，终点误差 %.3g m / %.3g s，T=%.9f s。\n', ...
    dObs,ePos,eT,pNom.tCross);

% ================= 嵌套回调 =================
    function updateFrame(tNow)
        tNow = min(max(tNow,0),tAll(end)); sNow = deval(pNom.sol,tNow)';
        set(hBall,'XData',sNow(1)+R*Xb,'YData',sNow(2)+R*Yb,'ZData',sNow(3)+R*Zb);
        set(hBallTop,'XData',sNow(1),'YData',sNow(2));
        % 过去实线：到 tNow 为止，末点追加当前点；未来虚线：从当前点到终点，不重复覆盖
        k = find(tAll <= tNow + 1e-12,1,'last'); if isempty(k), k = 1; end
        Spast = [[SAll(1:k,1);sNow(1)],[SAll(1:k,2);sNow(2)],[SAll(1:k,3);sNow(3)]];
        Sfut = [[sNow(1);SAll(k+1:end,1)],[sNow(2);SAll(k+1:end,2)],[sNow(3);SAll(k+1:end,3)]];
        set(hPlayed,'XData',Spast(:,1),'YData',Spast(:,2),'ZData',Spast(:,3));
        set(hFuture,'XData',Sfut(:,1),'YData',Sfut(:,2),'ZData',Sfut(:,3));
        set(hPlayedTop,'XData',Spast(:,1),'YData',Spast(:,2));
        set(hFutureTop,'XData',Sfut(:,1),'YData',Sfut(:,2));
        set(hSlider,'Value',tNow);
        [~,F] = ball_rhs(tNow,sNow(:),theta,cfg);
        set(hqG,'XData',sNow(1),'YData',sNow(2),'ZData',sNow(3),'UData',sc*F.gravity(1),'VData',sc*F.gravity(2),'WData',sc*F.gravity(3));
        set(hqD,'XData',sNow(1),'YData',sNow(2),'ZData',sNow(3),'UData',sc*F.drag(1),'VData',sc*F.drag(2),'WData',sc*F.drag(3));
        set(hqM,'XData',sNow(1),'YData',sNow(2),'ZData',sNow(3),'UData',sc*F.magnus(1),'VData',sc*F.magnus(2),'WData',sc*F.magnus(3));
        set(hTL,'String',sprintf(['t = %.3f s ｜ 位置 (%.2f, %.2f, %.2f) m ｜ ' ...
            '速度 (%.2f, %.2f, %.2f) m/s'],tNow,sNow(1:3),sNow(4:6)));
        set(hFL,'String',sprintf(['力 (N)：重力 (%+.2f, %+.2f, %+.2f) ｜ ' ...
            '阻力 (%+.2f, %+.2f, %+.2f) ｜ Magnus (%+.2f, %+.2f, %+.2f)'], ...
            F.gravity,F.drag,F.magnus));
    end

    function applyView(idx)
        % 只改变显示范围，不改变任何数据
        if idx == 1
            xr = [-0.4 11.6]; yr = [-0.4 3.8]; zr = [0 3];
        else
            xr = [7 11.6];    yr = [1.3 3.7]; zr = [1.2 2.8];
        end
        xlim(ax3,xr); ylim(ax3,yr); zlim(ax3,zr);
        pbaspect(ax3,[diff(xr) diff(yr) diff(zr)]);   % 1:1:1 单位比例，量程不被拉伸
        xlim(axTop,xr); ylim(axTop,yr);
        pbaspect(axTop,[diff(xr) diff(yr) 1]);        % 俯视等比例，量程保持不变
    end

    function cbView(src,~)
        applyView(get(src,'Value'));
    end

    function cbForce(src,~)
        onoff = 'off'; if get(src,'Value'), onoff = 'on'; end
        set([hqG hqD hqM],'Visible',onoff);
    end

    function cbPlay(~,~)
        if playing
            playing = false; stopOnly(timerH); set(hPlay,'String','播放'); return;
        end
        if get(hSlider,'Value') >= tAll(end)-1e-9, updateFrame(0); end
        t0 = get(hSlider,'Value'); tClock = tic; playing = true;
        start(timerH); set(hPlay,'String','暂停');
    end

    function cbTick(~,~)
        if ~playing, return; end
        tn = t0 + speed*toc(tClock);
        if tn >= tAll(end)
            updateFrame(tAll(end)); playing = false; stopOnly(timerH);
            set(hPlay,'String','重播'); return;
        end
        updateFrame(tn);
    end

    function cbReplay(~,~)
        if playing, stopOnly(timerH); playing = false; end
        t0 = 0; tClock = tic; updateFrame(0); playing = true;
        start(timerH); set(hPlay,'String','暂停');
    end

    function cbSlider(src,~)
        if playing, playing = false; stopOnly(timerH); set(hPlay,'String','播放'); end
        updateFrame(get(src,'Value'));
    end

    function cbSpeed(src,~)
        vv = [0.05 0.1 0.25 1]; speed = vv(get(src,'Value'));
        t0 = get(hSlider,'Value'); tClock = tic;
    end

    function cbCmp(src,~)
        onoff = 'off'; if get(src,'Value'), onoff = 'on'; end
        set([hG hD hGTop hDTop],'Visible',onoff);  % 三维与俯视图同步
    end

    function cbClose(~,~)
        stopTimer(timerH); delete(fig);
    end
end

% ================= 局部函数 =================
function p = predictPlane(theta,cfg,xPlane,Tmax)
%PREDICTPLANE ode45 0→Tmax；事件 x=xPlane(+1)、z=R(-1)，均终止
p = struct('sol',[],'te',[],'ye',[],'ie',[],'tCross',NaN,'status','timeout','y',NaN,'z',NaN);
opts = odeset('RelTol',cfg.relTol,'AbsTol',cfg.absTol,'MaxStep',cfg.odeMaxStep, ...
    'Events',@(t,s) planeEvents(t,s,cfg,xPlane));
p.sol = ode45(@(t,s) ball_rhs(t,s,theta,cfg),[0 Tmax],[cfg.r0;theta(1:3)'],opts);
p.te = p.sol.xe(:)'; p.ye = p.sol.ye'; p.ie = p.sol.ie(:)';
if isempty(p.te) || isempty(p.ie), return; end
p.tCross = p.te(1); p.status = 'ground';
if p.ie(1) == 1, p.status = 'plane'; p.y = p.ye(1,2); p.z = p.ye(1,3); end
end

function [value,isterminal,direction] = planeEvents(t,s,cfg,xPlane)
value = [s(1)-xPlane; s(3)-cfg.R];
if t < 1e-7, value(2) = abs(value(2)) + 1; end  % 起步初始接触不误终止
isterminal = [1; 1]; direction = [1; -1];
end

function T = tEndOf(p,Tmax)
if ~isempty(p.te) && isfinite(p.te(1)), T = p.te(1); else, T = Tmax; end
end

function [Fg,Fd,Fm] = forceSeries(ts,S,theta,cfg)
n = numel(ts); Fg = nan(n,3); Fd = Fg; Fm = Fg;
for k = 1:n
    [~,F] = ball_rhs(ts(k),S(k,:)',theta,cfg);
    Fg(k,:) = F.gravity'; Fd(k,:) = F.drag'; Fm(k,:) = F.magnus';
end
end

function [A,map] = captureRGB(fig,tmpDir,k)
try
    [A,map] = rgb2ind(frame2im(getframe(fig)),256);
catch  % getframe 不可靠时用临时 PNG，仅清理本函数生成的临时文件
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
