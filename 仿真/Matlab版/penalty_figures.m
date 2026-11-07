function penalty_figures()
%PENALTY_FIGURES 论文科学配图：系统流程、观测约束与轨迹分离、参数耦合与筛选敏感性。
% 仅读取 data/派生数据/罚球结果.json，不修改任何原始数据。
% 六张论文图按同一数据生成中英两版；白底，无图内标题。
% 英文在img/，中文在img/中文/；导出前后核对数据与坐标范围。

rootDir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
imgDir = fullfile(rootDir, 'img');
if ~exist(imgDir, 'dir'), mkdir(imgDir); end

familyFile = fullfile(rootDir, 'data', '派生数据', '罚球结果.json');
S = jsondecode(fileread(familyFile));
fam = S.trajectoryFamily;
mem = fam.members;
if iscell(mem), mem = [mem{:}]; end
ids = fam.exampleIndices(:)';
obsXY = S.observations.positions(:, 1:2);

% 只使用约定字段；初速由轨迹初速向量重建，避免依赖额外导出字段。
thetaAll = reshape([mem.theta], 6, []).';
speed = vecnorm(thetaAll(:, 1:3), 2, 2);
bz = thetaAll(:, 6);
wind = [mem.windY].';
density = [mem.density].';
Cd = [mem.Cd].';
dragProduct = density .* Cd;
accepted = logical([mem.accepted].');

blue = [0.05 0.32 0.65];
orange = [0.85 0.35 0.10];
darkGrey = [0.30 0.30 0.30];
lightGrey = [0.78 0.80 0.83];
fontName = 'Times New Roman';

pairReports={};
guards = cell(1, 6); % onCleanup 只关闭本函数自建图窗，不 close all

%% 图 1：系统与反演流程
fig1 = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 19.5 10.0]);
guards{1} = onCleanup(@() closeIfValid(fig1));
set(fig1, 'DefaultAxesFontName', fontName, 'DefaultAxesFontSize', 10, ...
    'DefaultTextFontName', fontName, 'DefaultTextFontSize', 10);
tl1 = tiledlayout(fig1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) 二维俯视系统示意
ax1 = nexttile(tl1);
hold(ax1, 'on');
plot(ax1, 11 * [1 1], [-3.66 3.66], 'k-', 'LineWidth', 2.0);
text(ax1, 11, 3.90, 'Goal plane x = 11 m', 'FontName', fontName, ...
    'FontSize', 9.5, 'HorizontalAlignment', 'center');
plot(ax1, 0, 0, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.15 0.15 0.15], ...
    'MarkerEdgeColor', 'k');
text(ax1, 0.10, -1.05, 'Penalty mark (0,0)', 'FontName', fontName, 'FontSize', 9.5);
plot(ax1, obsXY(:, 1), obsXY(:, 2), 'o', 'MarkerSize', 6.5, ...
    'MarkerFaceColor', orange, 'MarkerEdgeColor', 'k');
text(ax1, 0.15, 2.0, 'Observed centre t_1', ...
    'FontName', fontName, 'FontSize', 9.5);
text(ax1, 3.35, 1.0, 'Observed centre t_2', ...
    'FontName', fontName, 'FontSize', 9.5);
% 坐标箭头
quiver(ax1, 0, 0, 5.5, 0, 0, 'Color', darkGrey, 'LineWidth', 1.1, 'MaxHeadSize', 0.4);
quiver(ax1, 0, 0, 0, 2.7, 0, 'Color', darkGrey, 'LineWidth', 1.1, 'MaxHeadSize', 0.4);
text(ax1, 5.7, .12, 'x', 'FontName', fontName, 'FontSize', 10.5, 'Color', darkGrey);
text(ax1, 0.18, 2.85, 'y', 'FontName', fontName, 'FontSize', 10.5, 'Color', darkGrey);
% 风沿 +y
quiver(ax1, 6.0, 2.0, 0, 1.9, 0, 'Color', orange, 'LineWidth', 1.4, 'MaxHeadSize', 0.6);
text(ax1, 6.0, 4.10, 'Wind w_y', 'FontName', fontName, 'FontSize', 9.5, ...
    'Color', orange, 'HorizontalAlignment', 'center');
text(ax1, 12.55, -4.55, 'Schematic (top view, not to scale)', ...
    'FontName', fontName, 'FontSize', 8.5, 'FontAngle', 'italic', ...
    'HorizontalAlignment', 'right');
xlim(ax1, [-1.4 12.7]);
ylim(ax1, [-5.0 5.2]);
daspect(ax1, [1 1 1]);
axis(ax1, 'off');
labelPanel(ax1, '(a)', fontName);

% (b) 自上而下的反演流程框与箭头
ax2 = nexttile(tl1);
axis(ax2, 'off');
hold(ax2, 'on');
xlim(ax2, [0 1]);
ylim(ax2, [0 1]);
boxLabels = {sprintf('Given positions\n(t = 0.05 s, 0.10 s)'), ...
    sprintf('Prescribe\n(w_y, \\rho, C_D)'), ...
    sprintf('Fit\n(v_0, b_y, b_z)'), ...
    sprintf('Propagate to\ngoal / ground'), ...
    sprintf('Test residuals &\nphysical checks'), ...
    sprintf('Accepted trajectories')};
nBox = numel(boxLabels);
boxH = 0.120;
gap = 0.042;
boxW = 0.84;
boxX = 0.08;
for k = 1:nBox
    yTop = 1 - 0.03 - (k - 1) * (boxH + gap);
    yBot = yTop - boxH;
    if k == nBox
        faceColor = [0.90 0.94 0.99];
        edgeColor = blue;
    else
        faceColor = [0.96 0.97 0.99];
        edgeColor = [0.35 0.35 0.35];
    end
    rectangle(ax2, 'Position', [boxX yBot boxW boxH], 'Curvature', 0.16, ...
        'FaceColor', faceColor, 'EdgeColor', edgeColor, 'LineWidth', 0.9);
    text(ax2, boxX + boxW / 2, yBot + boxH / 2, boxLabels{k}, ...
        'FontName', fontName, 'FontSize', 10, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
    if k < nBox
        y1 = yBot - 0.009;
        y2 = yBot - gap + 0.006;
        quiver(ax2, 0.5, y1, 0, -(y1 - y2), 0, 'Color', greyArrow(), ...
            'LineWidth', 1.1, 'MaxHeadSize', 0.55, 'AutoScale', 'off');
    end
end
labelPanel(ax2, '(b)', fontName);

fig1File = fullfile(imgDir, '系统与反演流程.png');
pairReports{end+1}=exportPair(fig1,fig1File);
fprintf('导出完成：%s\n', fig1File);

%% 图 2：观测约束与轨迹分离
fig2 = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 19.5 8.4]);
guards{2} = onCleanup(@() closeIfValid(fig2));
set(fig2, 'DefaultAxesFontName', fontName, 'DefaultAxesFontSize', 10, ...
    'DefaultTextFontName', fontName, 'DefaultTextFontSize', 10);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

iA = ids(1);
iB = ids(2);

ax3 = nexttile(tl2);
hold(ax3, 'on');
fill(ax3, [0.5 6.5 6.5 0.5], [-5 -5 5 5], [0.88 0.92 0.97], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.55, 'HandleVisibility', 'off');
plot(ax3, [0.5 6.5], [5 5], '--', 'Color', [0.40 0.40 0.40], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax3, [0.5 6.5], [-5 -5], '--', 'Color', [0.40 0.40 0.40], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');
plot(ax3, [0.5 6.5], [0 0], '-', 'Color', [0.75 0.75 0.75], ...
    'LineWidth', 0.7, 'HandleVisibility', 'off');
plot(ax3, 1:6, rowResidualMm(mem(iA)), '-o', 'Color', blue, ...
    'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'MarkerSize', 5.5, ...
    'LineWidth', 1.3, 'DisplayName', sprintf('Representative A (%d)', iA));
plot(ax3, 1:6, rowResidualMm(mem(iB)), '-s', 'Color', orange, ...
    'MarkerFaceColor', orange, 'MarkerEdgeColor', 'k', 'MarkerSize', 5.5, ...
    'LineWidth', 1.3, 'DisplayName', sprintf('Representative B (%d)', iB));
xlim(ax3, [0.5 6.5]);
ylim(ax3, [-6 6]);
xticks(ax3, 1:6);
xticklabels(ax3, {'t_1 x', 't_1 y', 't_1 z', 't_2 x', 't_2 y', 't_2 z'});
xlabel(ax3, 'Observation coordinate');
ylabel(ax3, 'Residual (mm)');
styleAxes(ax3, fontName);
legend(ax3, 'Location', 'southoutside', 'FontName', fontName, 'FontSize', 9, 'Box', 'off');
labelPanel(ax3, '(a)', fontName);

ax4 = nexttile(tl2);
hold(ax4, 'on');
hFam = [];
for j = find(accepted).'
    r = mem(j).trajectory;
    if isempty(r), continue; end
    r = r(r(:, 1) >= 8.5, :);
    if size(r, 1) < 2, continue; end
    if isempty(hFam)
        hFam = plot(ax4, r(:, 1), r(:, 2), '-', 'Color', lightGrey, ...
            'LineWidth', 0.7, 'DisplayName', 'Accepted trajectories');
    else
        plot(ax4, r(:, 1), r(:, 2), '-', 'Color', lightGrey, ...
            'LineWidth', 0.7, 'HandleVisibility', 'off');
    end
end
tailA = mem(iA).trajectory;
tailA = tailA(tailA(:, 1) >= 8.5, :);
tailB = mem(iB).trajectory;
tailB = tailB(tailB(:, 1) >= 8.5, :);
hA = plot(ax4, tailA(:, 1), tailA(:, 2), '-', 'Color', blue, ...
    'LineWidth', 1.8, 'DisplayName', 'Representative A');
hB = plot(ax4, tailB(:, 1), tailB(:, 2), '-', 'Color', orange, ...
    'LineWidth', 1.8, 'DisplayName', 'Representative B');
xline(ax4, 11, ':', 'Color', [0.55 0.55 0.55], 'HandleVisibility', 'off');
scatter(ax4, 11, mem(iA).endpoint(2), 55, blue, 'filled', ...
    'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
scatter(ax4, 11, mem(iB).endpoint(2), 55, orange, 'filled', ...
    'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
text(ax4, 10.90, mem(iA).endpoint(2), 'A', 'FontName', fontName, ...
    'FontSize', 10.5, 'FontWeight', 'bold', 'Color', blue, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
text(ax4, 10.90, mem(iB).endpoint(2), 'B', 'FontName', fontName, ...
    'FontSize', 10.5, 'FontWeight', 'bold', 'Color', orange, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
yAll = [tailA(:, 2); tailB(:, 2); mem(iA).endpoint(2); mem(iB).endpoint(2)];
for j = find(accepted).'
    r = mem(j).trajectory;
    if ~isempty(r)
        r = r(r(:, 1) >= 8.5, :);
        yAll = [yAll; r(:, 2)]; %#ok<AGROW>
    end
end
yPad = max(0.08 * (max(yAll) - min(yAll)), 0.02);
ylim(ax4, [min(yAll) - yPad, max(yAll) + yPad]);
xlim(ax4, [8.5 11.15]);
xlabel(ax4, 'x (m)');
ylabel(ax4, 'y (m)');
styleAxes(ax4, fontName);
if isempty(hFam)
    legend(ax4, [hA, hB], {'Representative A', 'Representative B'}, ...
        'Location', 'southoutside', 'FontName', fontName, 'FontSize', 9, 'Box', 'off');
else
    legend(ax4, [hFam, hA, hB], ...
        {'Accepted trajectories', 'Representative A', 'Representative B'}, ...
        'Location', 'southoutside', 'FontName', fontName, 'FontSize', 9, 'Box', 'off');
end
labelPanel(ax4, '(b)', fontName);

fig2File = fullfile(imgDir, '观测约束与轨迹分离.png');
pairReports{end+1}=exportPair(fig2,fig2File);
fprintf('导出完成：%s\n', fig2File);

%% 图 3：参数耦合与筛选敏感性
fig3 = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 19.8 8.0]);
guards{3} = onCleanup(@() closeIfValid(fig3));
set(fig3, 'DefaultAxesFontName', fontName, 'DefaultAxesFontSize', 10, ...
    'DefaultTextFontName', fontName, 'DefaultTextFontSize', 10);
tl3 = tiledlayout(fig3, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) 三个风切片内的有效阻力与拟合初速
ax5 = nexttile(tl3);
hold(ax5, 'on');
wSlices = [-8 0 8];
wColors = [blue; orange; darkGrey];
hWind = gobjects(1, numel(wSlices));
for k = 1:numel(wSlices)
    okSlice = accepted & abs(wind - wSlices(k)) < 1e-9;
    if any(okSlice)
        hWind(k) = scatter(ax5, dragProduct(okSlice), speed(okSlice), 28, ...
            wColors(k, :), 'filled', 'MarkerEdgeColor', [0.15 0.15 0.15], ...
            'LineWidth', 0.4, 'DisplayName', sprintf('w_y = %+g m/s', wSlices(k)));
    else
        hWind(k) = scatter(ax5, NaN, NaN, 28, wColors(k, :), 'filled', ...
            'DisplayName', sprintf('w_y = %+g m/s', wSlices(k)));
    end
end
rej = ~accepted & any(abs(wind - wSlices) < 1e-9, 2);
hRej = scatter(ax5, dragProduct(rej), speed(rej), 24, [0.62 0.62 0.62], ...
    'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.45 0.45 0.45], ...
    'LineWidth', 0.7, 'DisplayName', 'Rejected (residual)');
xlabel(ax5, 'Effective drag \rho C_D (kg m^{-3})');
ylabel(ax5, 'Fitted initial speed (m s^{-1})');
styleAxes(ax5, fontName);
legend(ax5, [hWind, hRej], 'Location', 'southoutside', 'FontName', fontName, ...
    'FontSize', 8.5, 'Box', 'off');
xlim(ax5,[min(dragProduct)-.03,max(dragProduct)+.03]);ylim(ax5,[min(speed)-.04,max(speed)+.04]);
labelPanel(ax5, '(a)', fontName);

% (b) 所有合格点的指定风与拟合 b_z
ax6 = nexttile(tl3);
hold(ax6, 'on');
scatter(ax6, wind(accepted), bz(accepted), 28, blue, 'filled', ...
    'MarkerEdgeColor', [0.10 0.10 0.10], 'LineWidth', 0.5);
xlabel(ax6, 'Prescribed wind w_y (m s^{-1})');
ylabel(ax6, {'Fitted b_z (s^{-1})'; '(effective, not measured spin)'});
styleAxes(ax6, fontName);
if any(accepted)
    xlim(ax6, [min(wind(accepted)) - 1, max(wind(accepted)) + 1]);
    yb = bz(accepted);
    yPad = max(0.08 * (max(yb) - min(yb)), 0.02);
    ylim(ax6, [min(yb) - yPad, max(yb) + yPad]);
end
labelPanel(ax6, '(b)', fontName);set(findobj(ax6,'Type','text','String','(b)'),'Position',[.8 .98 0]);

% (c) 残差容差与接受数（计数，不是概率）
ax7 = nexttile(tl3);
toleranceMm = [fam.summary.tolerance] * 1000;
acceptedCounts = [fam.summary.nAccepted];
plot(ax7, toleranceMm, acceptedCounts, '-o', 'Color', blue, ...
    'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
    'LineWidth', 1.6);
xlim(ax7, [0, max(toleranceMm) + 1]);
ylim(ax7, [0, 1.06*fam.nCandidates]);
xticks(ax7, toleranceMm);
xlabel(ax7, 'Residual tolerance (mm)');
ylabel(ax7, 'Accepted candidates (count)');
styleAxes(ax7, fontName);
labelPanel(ax7, '(c)', fontName);

fig3File = fullfile(imgDir, '参数耦合与筛选敏感性.png');
pairReports{end+1}=exportPair(fig3,fig3File);

fprintf('导出完成：%s\n', fig3File);

% 与正文/附录共用同一数据的三张统计图。
fig4=figure('Visible','off','Color','w','Units','centimeters','Position',[2 2 16 11]);
guards{4}=onCleanup(@()closeIfValid(fig4));ax=axes(fig4);hold(ax,'on');
points=reshape([mem(accepted).endpoint],3,[]);
scatter(ax,points(2,:),points(3,:),22,[mem(accepted).Cd],'filled');
cb=colorbar(ax);cb.Label.String='Drag coefficient C_D';
xlabel(ax,'y at x = 11 m (m)');ylabel(ax,'z at x = 11 m (m)');
styleAxes(ax,fontName);axis(ax,'equal');xlim(ax,[min(points(2,:))-.1,max(points(2,:))+.1]);ylim(ax,[min(points(3,:))-.1,max(points(3,:))+.1]);
pairReports{end+1}=exportPair(fig4,fullfile(imgDir,'全程约束球门散点.png'));
fid=fopen(fullfile(rootDir,'data','派生数据','论文配图核查.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(struct('pairs',{pairReports},'sameDataAndAxes',true),'PrettyPrint',true));fclose(fid);
end

% -------------------------------------------------------------------------
function r = rowResidualMm(s)
% 残差为 2x3（两个观测时刻 x/y/z），逐行展开为 6 个坐标残差，单位 mm。
res = s.residual;
r = [res(1, :), res(2, :)] * 1000;
end

function styleAxes(ax, fontName)
set(ax, 'FontName', fontName, 'FontSize', 10, 'Box', 'on', ...
    'LineWidth', 0.6, 'TickDir', 'out');
grid(ax, 'on');
end

function labelPanel(ax, txt, fontName)
text(ax, 0.01, .98, txt, 'Units', 'normalized', 'FontName', fontName, ...
    'FontSize', 11, 'FontWeight', 'bold', 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left');
end

function c = greyArrow()
c = [0.35 0.35 0.35];
end

function closeIfValid(fig)
if isgraphics(fig)
    close(fig);
end
end


function report=exportPair(fig,path)
% 同一图形先导出英文，再只替换文字导出中文；数据与坐标范围逐项核对。
geometry=graphGeometry(fig);
exportgraphics(fig,path,'Resolution',220);
objects=findall(fig);
for k=1:numel(objects)
 h=objects(k);
 if isprop(h,'FontName'),h.FontName='SimSun';end
 if strcmp(h.Type,'colorbar'),h.Label.String=translateLabel(h.Label.String);end
 if isprop(h,'String') && any(strcmp(h.Type,{'text','legend'}))
  h.String=translateLabel(h.String);
 end
end
assert(isequaln(geometry,graphGeometry(fig)),'双语转换不得改变数据或坐标范围。');
[parent,name,ext]=fileparts(path);cn=fullfile(parent,'中文');if ~exist(cn,'dir'),mkdir(cn);end
cnpath=fullfile(cn,[name ext]);exportgraphics(fig,cnpath,'Resolution',220);
report=struct('english',path,'chinese',cnpath,'sameDataAndAxes',true);
end
function g=graphGeometry(fig)
g={};objects=findall(fig);
for k=1:numel(objects)
 h=objects(k);
 for field={'XData','YData','ZData','XLim','YLim','ZLim'}
  if isprop(h,field{1}),g{end+1}=h.(field{1});end %#ok<AGROW>
 end
end
end
function out=translateLabel(in)
if iscell(in),out=cellfun(@translateLabel,in,'UniformOutput',false);return;end
if isstring(in),out=arrayfun(@translateLabel,in,'UniformOutput',false);return;end
out=in;
pairs={...
 'Goal plane x = 11 m','球门平面 x = 11 m';...
 'Penalty mark (0,0)','罚球点 (0,0)';...
 'Observed centre','观测球心';...
 'Wind w_y','横风 w_y';...
 'Schematic (top view, not to scale)','水平几何示意（非比例图）';...
 'Given positions','题给位置';...
 'Prescribe','指定环境';...
 'Fit','反求';...
 'Propagate to','全程积分至';...
 'goal / ground','球门 / 触地';...
 'Test residuals &','检验观测残差及';...
 'physical checks','物理约束';...
 'Accepted trajectories','符合观测的轨迹';...
 'Representative A','代表解甲';...
 'Representative B','代表解乙';...
 'Observation coordinate','观测坐标';...
 'Residual (mm)','坐标残差 (mm)';...
 'Rejected (residual)','残差未通过';...
 'Effective drag','有效阻力参数';...
 'Fitted initial speed','反求初速度';...
 'Prescribed wind','指定横风';...
 'Fitted b_z','反求 b_z';...
 '(effective, not measured spin)','（有效系数）';...
 'Residual tolerance (mm)','观测容差 (mm)';...
 'Accepted candidates (count)','符合观测的轨迹数（个）';...
 'Drag coefficient C_D','阻力系数 C_D';...
 'y at x = 11 m (m)','球门平面 y (m)';...
 'z at x = 11 m (m)','球门平面 z (m)';...
 'Bootstrap predictions','自助重拟合落点';...
 '95% joint ellipse','95% 联合椭圆';...
 'Conditional point','条件参考点';...
 'Goal frame','球门框';...
 'Reference','参考点';...
 'Change in y (cm)','y 方向偏移 (cm)';...
 'Change in z (cm)','z 方向偏移 (cm)';...
 'wind_y=','横风=';...
 'Cd exponent p=','阻力速度指数 p=';...
 'b decay tau=','侧力衰减时间='};
% 长标签先替换，避免Fit误替换Fitted。
[~,ix]=sort(cellfun(@numel,pairs(:,1)),'descend');
for k=ix(:)',out=strrep(out,pairs{k,1},pairs{k,2});end
if strcmp(out,'A'),out='甲';elseif strcmp(out,'B'),out='乙';end
end
