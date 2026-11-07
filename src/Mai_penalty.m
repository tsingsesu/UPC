function result = Mai_penalty(mode, obs)
%MAI_PENALTY  FIFA点球早期观测反演、球门落点预测与不确定性评估
%
%   result = Mai_penalty(mode)
%
%   mode = 'full'  （默认）完整Bootstrap与合成反演样本
%   mode = 'quick' （少量样本，用于快速流程检查）
%
%   物理状态 r'=v，v' = g - K*Cd*||u||*u + cross(b,u)，u = v - w。
%   b = [bx;by;bz] 是有效Magnus系数向量，单位 1/s；bx固定为0，
%   因此本程序不声称能够唯一恢复真实自旋。
%   theta = [vx vy vz Cd by bz]。
%
%   本文件只做实现，不在此处运行；运行与结果核验由Codex完成。
%   所有正式派生数据写入项目根下 data/，图片写入 img/。

if nargin < 1 || isempty(mode)
    mode = 'full';
end
mode = lower(char(mode));
if ~any(strcmp(mode, {'full','quick'}))
    error('mode 必须是 ''full'' 或 ''quick''。');
end
tStart = tic;

seed = 20261107;
rng(seed, 'twister');

% ---------- 项目路径（由当前文件位置计算，不含个人绝对路径） ----------
thisFile = mfilename('fullpath');
rootDir = fileparts(fileparts(thisFile));
dataDir = fullfile(rootDir, 'data', '派生数据');
imgDir  = fullfile(rootDir, 'img');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
if ~exist(imgDir,  'dir'), mkdir(imgDir);  end

% ---------- 固定物理与数值参数 ----------
cfg = struct();
cfg.R       = 0.11;                 % 球半径 m
cfg.m       = 0.43;                 % 球质量 kg
cfg.rho     = 1.225;                % 空气密度 kg/m^3
cfg.g       = [0; 0; -9.80665];     % 重力加速度 m/s^2
cfg.r0      = [0; 0; cfg.R];        % 发球瞬间球心位置 m
cfg.K       = cfg.rho*pi*cfg.R^2/(2*cfg.m);
cfg.w       = [0; 0; 0];            % 风速 m/s
cfg.bx      = 0;                    % 有效Magnus系数x分量固定为0
cfg.h       = 0.005;                % 观测拟合固定RK4步长 s
cfg.cdExponent = 0;                 % Cd=Cd0*(speed/30)^p，p=0即常数
cfg.cdRefSpeed = 30;                % 变Cd参考速度 m/s
cfg.bDecayTau  = inf;               % b(t)=b0*exp(-t/tau)，inf即不衰减
cfg.odeMaxStep = 0.01;              % ode45预测最大步长 s
cfg.sigma   = 0.01;                 % 观测各坐标独立高斯噪声标准差 m
cfg.lb      = [10, -20, -5, 0.05, -1.5, -1.5];
cfg.ub      = [50,  20, 25, 0.60,  1.5,  1.5];
cfg.freeMask   = true(1, 6);
cfg.thetaFixed = zeros(1, 6);

obsT = [0.05, 0.10];
if nargin<2 || isempty(obs)
    obs=[1.50,.45,.44;2.96,.86,.74];
end
validateattributes(obs,{'numeric'},{'size',[2,3],'finite','real'});

fprintf('[Mai_penalty] 模式=%s，开始完整反演。\n', mode);

% ---------- 1. 完整有效Magnus模型拟合（多初值核查局部解） ----------
starts = [24, 3, 5, 0.25,  0.0,  0.0;
          30, 6, 8, 0.30,  0.0, -0.3;
          20, 0, 3, 0.20,  0.2,  0.2;
          40,10,12, 0.40, -0.2,  0.4];
[fitBest, fitStarts] = fitWithStarts(starts, obsT, obs, cfg);
if ~isfinite(fitBest.sse) || fitBest.exitflag<=0
    error('完整模型的所有初值均未得到有效拟合，请检查输入与工具箱。');
end
thetaHat = fitBest.theta;
fprintf('[Mai_penalty] 拟合完成：RMSE=%.6g m，退出标志=%d，边界命中=%d。\n', ...
    fitBest.rmse, fitBest.exitflag, fitBest.boundaryHit);

% 基准预测（x=11 球门平面）与球心全越线 x=11.11 的几何敏感性
predBase   = predictGoal(thetaHat, cfg, 11.00);
predCenter = predictGoal(thetaHat, cfg, 11.11);
fprintf('[Mai_penalty] 基准预测状态=%s，y=%.6g，z=%.6g。\n', ...
    predBase.status, predBase.yGoal, predBase.zGoal);

% ---------- 2. 观测位置雅可比、参数协方差与预测协方差 ----------
[J, Jinfo] = obsJacobian(thetaHat, cfg, obsT);
[~, S, ~] = svd(J, 'econ');
sv = diag(S);
condJ = cond(J);
covTheta = cfg.sigma^2 * pinv(J'*J);
[G, Ginfo] = goalJacobian(thetaHat, cfg);
if all(isfinite(G(:)))
    C = G * covTheta * G';
else
    C = nan(2, 2);
end
ell = ellipseFromCov(C);
marginalHalfWidth = 1.959964 * sqrt(max(diag(C), 0));

unc = struct();
unc.J = J;
unc.J_perturbation = Jinfo.perturbation;
unc.singularValues = sv;
unc.conditionNumber = condJ;
unc.covTheta = covTheta;
unc.G = G;
unc.G_perturbation = Ginfo.perturbation;
unc.C = C;
unc.joint95Axes = ell.axes;
unc.joint95Vectors = ell.vectors;
unc.marginal95HalfWidth_yz = marginalHalfWidth;
unc.note = ['参数协方差按 sigma^2*pinv(J''*J) 计算，J 是观测位置对参数的雅可比；' ...
    '预测映射 g(theta)=[yGoal,zGoal] 的数值雅可比为 G，C=G*covTheta*G''。' ...
    '区间仅为局部高斯线性化结果，不处理参数边界，也不代表无条件概率众数。'];

% ---------- 3. 不同观测噪声水平情景（不重新拟合，仅传播噪声假设） ----------
sigmas = [0.002, 0.005, 0.010, 0.020];
noiseScenarios = repmat(struct('sigma',NaN, 'covTheta',[], 'C',[], ...
    'joint95Axes',[], 'marginal95HalfWidth_yz',[]), numel(sigmas), 1);
for k = 1:numel(sigmas)
    covS = sigmas(k)^2 * pinv(J'*J);
    CS = G * covS * G';
    eS = ellipseFromCov(CS);
    noiseScenarios(k).sigma = sigmas(k);
    noiseScenarios(k).covTheta = covS;
    noiseScenarios(k).C = CS;
    noiseScenarios(k).joint95Axes = eS.axes;
    noiseScenarios(k).marginal95HalfWidth_yz = 1.959964 * sqrt(max(diag(CS), 0));
end
fprintf('[Mai_penalty] 噪声情景完成：%d 个假设噪声水平。\n', numel(sigmas));

% ---------- 4. 四个基线模型 ----------
baselines = repmat(struct('name','', 'theta',nan(1,6), 'v0',nan(1,3), ...
    'accel',nan(1,3), 'rmse',NaN, 'residual',nan(1,6), 'status','', ...
    'goal',false, 'yGoal',NaN, 'zGoal',NaN, 'tCross',NaN), 4, 1);

% 基线1：仅重力弹道，v0由两次观测线性最小二乘
v0G = zeros(1, 3);
for c = 1:3
    v0G(c) = sum(obsT .* (obs(:,c)' - cfg.r0(c) - 0.5*cfg.g(c)*obsT.^2)) / sum(obsT.^2);
end
aG = cfg.g';
P1 = repmat(cfg.r0', 2, 1) + obsT(:)*v0G + 0.5*(obsT(:).^2)*aG;
p1 = analyticPrediction(cfg.r0', v0G, aG, cfg);
baselines(1) = makeBaseline('重力弹道（v0线性最小二乘）', nan(1,6), v0G, aG, P1, obs, p1);

% 基线2：常加速度恰好插值两次观测，解析求x=11最早正实根
v0C = (4*obs(1,:) - obs(2,:) - 3*cfg.r0') / 0.10;
aC  = (obs(2,:) - 2*obs(1,:) + cfg.r0') / 0.05^2;
P2 = repmat(cfg.r0', 2, 1) + obsT(:)*v0C + 0.5*(obsT(:).^2)*aC;
p2 = analyticPrediction(cfg.r0', v0C, aC, cfg);
baselines(2) = makeBaseline('常加速度精确插值', nan(1,6), v0C, aC, P2, obs, p2);

% 基线3：仅阻力Cd与v0
cfg3 = cfg;
cfg3.freeMask = [true, true, true, true, false, false];
cfg3.thetaFixed = [0, 0, 0, 0, 0, 0];
[fit3, ~] = fitWithStarts([24, 3, 5, 0.25, 0, 0], obsT, obs, cfg3);
p3 = predictGoal(fit3.theta, cfg, 11.00);
baselines(3) = makeBaseline('仅阻力（Cd与v0）', fit3.theta, nan(1,3), nan(1,3), ...
    predictPositions(fit3.theta, obsT, cfg), obs, p3);
baselines(3).rmse = fit3.rmse;
baselines(3).residual = fit3.resid;

% 基线4：完整有效Magnus模型
baselines(4) = makeBaseline('完整有效Magnus模型', thetaHat, nan(1,3), nan(1,3), ...
    predictPositions(thetaHat, obsT, cfg), obs, predBase);
baselines(4).rmse = fitBest.rmse;
baselines(4).residual = fitBest.resid;
fprintf('[Mai_penalty] 四个基线模型完成。\n');

% ---------- 5. 参数Bootstrap（重新拟合） ----------
if strcmp(mode, 'full'), nb = 600; else, nb = 40; end
rng(seed, 'twister');
bsSamples = repmat(struct('index',0, 'theta',nan(1,6), 'ok',false, ...
    'exitflag',0, 'boundaryHit',false, 'status','', 'goal',false, ...
    'yGoal',NaN, 'zGoal',NaN, 'sse',NaN), nb, 1);
for b = 1:nb
    obsB = predictPositions(thetaHat, obsT, cfg) + cfg.sigma * randn(size(obs));
    fitB = fitWithStarts(thetaHat, obsT, obsB, cfg);
    bsSamples(b).index = b;
    bsSamples(b).theta = fitB.theta;
    bsSamples(b).exitflag = fitB.exitflag;
    bsSamples(b).boundaryHit = fitB.boundaryHit;
    bsSamples(b).sse = fitB.sse;
    okB = fitB.exitflag > 0 && all(isfinite(fitB.theta));
    bsSamples(b).ok = okB;
    if okB
        predB = predictGoal(fitB.theta, cfg, 11.00);
        bsSamples(b).status = predB.status;
        bsSamples(b).goal = predB.goal;
        bsSamples(b).yGoal = predB.yGoal;
        bsSamples(b).zGoal = predB.zGoal;
    else
        bsSamples(b).status = 'fit_fail';
    end
    if mod(b, 50) == 0 || b == nb
        fprintf('[Mai_penalty] Bootstrap 进度 %d/%d。\n', b, nb);
    end
end

okBS = [bsSamples.ok]';
thBS = vertcat(bsSamples.theta);
validBS = okBS & all(isfinite(thBS), 2);
if any(validBS)
    qBS = prctile(thBS(validBS,:), [2.5, 50, 97.5], 1);
    meanBS = mean(thBS(validBS,:), 1);
    medBS = median(thBS(validBS,:), 1);
    stdBS = std(thBS(validBS,:), 0, 1);
else
    qBS = nan(3,6); meanBS = nan(1,6); medBS = nan(1,6); stdBS = nan(1,6);
end
goalBS = [bsSamples.goal]' & okBS;
bsSummary = struct();
bsSummary.nSamples = nb;
bsSummary.nSuccessful = sum(okBS);
bsSummary.nFailed = sum(~okBS);
bsSummary.nBoundary = sum([bsSamples.boundaryHit]');
bsSummary.nValidForQuantiles = sum(validBS);
bsSummary.goalStabilityFrequency = sum(goalBS) / max(1, sum(okBS));
bsSummary.goalFractionAllDraws = sum(goalBS) / nb;
bsSummary.quantiles_2_5_50_97_5 = qBS;
bsSummary.mean = meanBS;
bsSummary.median = medBS;
bsSummary.std = stdBS;
bsSummary.parameterOrder = {'vx','vy','vz','Cd','by','bz'};
bsSummary.note = ['goalStabilityFrequency 是Bootstrap重拟合中预测越过球门的频率，' ...
    '不是实际进球概率；失败样本不计入其分母，另给全部抽样分母的比例。' ...
    '不能把筛选后的进球样本区间冒充无条件区间。'];
planeMask=okBS & isfinite([bsSamples.yGoal]') & isfinite([bsSamples.zGoal]');
bsSummary.nPlaneCrossings=sum(planeMask);
bsSummary.nGround=sum(strcmp({bsSamples.status},'ground'));
bsSummary.nOutside=sum(strcmp({bsSamples.status},'cross_out'));
bsSummary.goalPointQuantiles=prctile([[bsSamples(planeMask).yGoal]',[bsSamples(planeMask).zGoal]'],[2.5,50,97.5],1);
bsSummary.goalPointMean=mean([[bsSamples(planeMask).yGoal]',[bsSamples(planeMask).zGoal]'],1);
bsSummary.goalPointStd=std([[bsSamples(planeMask).yGoal]',[bsSamples(planeMask).zGoal]'],0,1);
fprintf('[Mai_penalty] Bootstrap完成：成功 %d/%d，进球稳定性频率 %.3f。\n', ...
    bsSummary.nSuccessful, nb, bsSummary.goalStabilityFrequency);

% ---------- 6. 数值验证 ----------
ver = struct();

% 6.1 RK4在Cd=b=0时与解析重力弹道比较
cfgV = cfg;
thetaV = [25, 4, 6, 0, 0, 0];
tV = 0.05:0.005:0.20;
PV = predictPositions(thetaV, tV, cfgV);
PA = repmat(cfg.r0', numel(tV), 1) + tV(:)*thetaV(1:3) + ...
    0.5*(tV(:).^2)*cfg.g';
ver.rk4GravityMaxError = max(abs(PV(:) - PA(:)));

% 6.2 无重力、无Magnus时一维阻力解析解
cfgD = cfg;
cfgD.g = [0;0;0];
thetaD = [18, 0, 0, 0.25, 0, 0];
tD = 0.05:0.005:0.20;
PD = predictPositions(thetaD, tD, cfgD);
xD_analytic = log(1 + cfg.K*thetaD(4)*thetaD(1)*tD) / (cfg.K*thetaD(4));
ver.dragAnalyticMaxError = max(abs(PD(:,1)' - xD_analytic));

% 6.3 步长与ode45最大步长比较
hs = [0.010, 0.005, 0.0025];
ver.stepSizes = repmat(struct('h',NaN, 'obsRmse',NaN, 'status','', ...
    'yGoal',NaN, 'zGoal',NaN), numel(hs), 1);
for k = 1:numel(hs)
    cfgH = cfg; cfgH.h = hs(k);
    PH = predictPositions(thetaHat, obsT, cfgH);
    predH = predictGoal(thetaHat, cfgH, 11.00);
    ver.stepSizes(k).h = hs(k);
    ver.stepSizes(k).obsRmse = sqrt(mean((PH(:)-obs(:)).^2));
    ver.stepSizes(k).status = predH.status;
    ver.stepSizes(k).yGoal = predH.yGoal;
    ver.stepSizes(k).zGoal = predH.zGoal;
end
ms = [0.020, 0.010, 0.005];
ver.odeMaxSteps = repmat(struct('maxStep',NaN, 'status','', ...
    'yGoal',NaN, 'zGoal',NaN, 'diffFromFinest',NaN), numel(ms), 1);
for k = 1:numel(ms)
    cfgM = cfg; cfgM.odeMaxStep = ms(k);
    predM = predictGoal(thetaHat, cfgM, 11.00);
    ver.odeMaxSteps(k).maxStep = ms(k);
    ver.odeMaxSteps(k).status = predM.status;
    ver.odeMaxSteps(k).yGoal = predM.yGoal;
    ver.odeMaxSteps(k).zGoal = predM.zGoal;
end
ref = [ver.odeMaxSteps(end).yGoal; ver.odeMaxSteps(end).zGoal];
for k = 1:numel(ms)
    if all(isfinite([ver.odeMaxSteps(k).yGoal, ver.odeMaxSteps(k).zGoal])) && all(isfinite(ref))
        ver.odeMaxSteps(k).diffFromFinest = norm([ver.odeMaxSteps(k).yGoal; ...
            ver.odeMaxSteps(k).zGoal] - ref);
    end
end

% 6.4 无风Magnus功率与阻力耗能率检查
tPow = 0.10:0.01:0.50;
SPow = simulateAt(thetaHat, tPow, cfg, 0.01);
bPow = [cfg.bx; thetaHat(5); thetaHat(6)];
magnusDot = zeros(numel(tPow),1);
dragPower = zeros(numel(tPow),1);
for k = 1:numel(tPow)
    v = SPow(k,4:6)';
    magnusDot(k) = v' * cross(bPow, v);
    dragPower(k) = -cfg.K*thetaHat(4)*norm(v)*(v'*v);
end
ver.magnusPowerMaxAbs = max(abs(magnusDot));
ver.dragPowerMax = max(dragPower);
ver.magnusPowerPass = ver.magnusPowerMaxAbs < 1e-10;
ver.dragPowerPass = ver.dragPowerMax <= 1e-12;

% 6.5 独立合成轨迹反演检验
if strcmp(mode, 'full'), nSyn = 150; else, nSyn = 10; end
rng(seed, 'twister');
startsSyn = [28, 6, 6, 0.25,  0.0,  0.0;
             24, 4, 5, 0.20,  0.1, -0.2];
synSamples = repmat(struct('index',0, 'thetaTrue',nan(1,6), 'thetaHat',nan(1,6), ...
    'ok',false, 'exitflag',0, 'boundaryHit',false, 'fitStatus','', ...
    'truthStatus','', 'truthGoal',false, 'fitGoal',false, ...
    'yTrue',NaN, 'zTrue',NaN, 'yHat',NaN, 'zHat',NaN, ...
    'coverValid',false, 'covered',false, 'covPred',nan(2,2), ...
    'predSqErr',NaN), nSyn, 1);
coverThreshold = chi2inv(0.95, 2);
for k = 1:nSyn
    thetaTrue = [24 + 10*rand, 3 + 7*rand, 4 + 5*rand, ...
        0.15 + 0.25*rand, -0.15 + 0.30*rand, -0.60 + 0.70*rand];
    STrue = simulateAt(thetaTrue, obsT, cfg, 1e-3);
    obsSyn = STrue(:,1:3) + cfg.sigma*randn(2,3);
    fitS = fitWithStarts(startsSyn, obsT, obsSyn, cfg);
    cfgTruth=cfg; cfgTruth.h=.0005; cfgTruth.odeMaxStep=.001;
    predTrue = predictGoal(thetaTrue, cfgTruth, 11.00);
    predHat  = predictGoal(fitS.theta, cfg, 11.00);
    synSamples(k).index = k;
    synSamples(k).thetaTrue = thetaTrue;
    synSamples(k).thetaHat = fitS.theta;
    synSamples(k).exitflag = fitS.exitflag;
    synSamples(k).boundaryHit = fitS.boundaryHit;
    synSamples(k).truthStatus = predTrue.status;
    synSamples(k).truthGoal = predTrue.goal;
    synSamples(k).yTrue = predTrue.yGoal;
    synSamples(k).zTrue = predTrue.zGoal;
    okS = fitS.exitflag > 0 && all(isfinite(fitS.theta));
    synSamples(k).ok = okS;
    if okS
        synSamples(k).fitStatus = predHat.status;
        synSamples(k).fitGoal = predHat.goal;
        synSamples(k).yHat = predHat.yGoal;
        synSamples(k).zHat = predHat.zGoal;
        if all(isfinite([predTrue.yGoal,predTrue.zGoal,predHat.yGoal,predHat.zGoal]))
            synSamples(k).predSqErr = norm([predHat.yGoal-predTrue.yGoal; ...
                predHat.zGoal-predTrue.zGoal])^2;
        end
        Js = obsJacobian(fitS.theta, cfg, obsT);
        Gs = goalJacobian(fitS.theta, cfg);
        if all(isfinite(Gs(:)))
            covP = Gs * (cfg.sigma^2*pinv(Js'*Js)) * Gs';
            synSamples(k).covPred = covP;
            if all(isfinite([predTrue.yGoal,predTrue.zGoal,predHat.yGoal,predHat.zGoal])) && all(isfinite(covP(:)))
                delta = [predHat.yGoal - predTrue.yGoal; ...
                         predHat.zGoal - predTrue.zGoal];
                d2 = delta' * pinv(covP) * delta;
                synSamples(k).coverValid = true;
                synSamples(k).covered = d2 <= coverThreshold;
            end
        end
    else
        synSamples(k).fitStatus = 'fit_fail';
    end
    if mod(k, 25) == 0 || k == nSyn
        fprintf('[Mai_penalty] 合成反演进度 %d/%d。\n', k, nSyn);
    end
end

okSyn = [synSamples.ok]';
truthGoalSyn = [synSamples.truthGoal]';
coverValidSyn = [synSamples.coverValid]';
coveredSyn = [synSamples.covered]';
thTrue = vertcat(synSamples.thetaTrue);
thHatS = vertcat(synSamples.thetaHat);
validSyn = okSyn & all(isfinite(thHatS), 2);
if any(validSyn)
    paramRmse = sqrt(mean((thHatS(validSyn,:) - thTrue(validSyn,:)).^2, 1));
else
    paramRmse = nan(1,6);
end
predRmse = nan;
tmpPred = [synSamples.predSqErr]';
if any(isfinite(tmpPred))
    predRmse = sqrt(mean(tmpPred(isfinite(tmpPred))));
end
synSummary = struct();
synSummary.nSamples = nSyn;
synSummary.nSuccessful = sum(okSyn);
synSummary.nFailed = sum(~okSyn);
synSummary.nBoundary = sum([synSamples.boundaryHit]');
synSummary.parameterRmse = paramRmse;
synSummary.parameterOrder = {'vx','vy','vz','Cd','by','bz'};
synSummary.predictionRmse_planeSamples = predRmse;
synSummary.coverageDenominator = sum(coverValidSyn);
synSummary.coverageRate = sum(coveredSyn & coverValidSyn) / max(1, sum(coverValidSyn));
synSummary.coverageThresholdChi2_95_2 = coverThreshold;
synSummary.nTruthPlaneCrossings=sum(isfinite([synSamples.yTrue]));
synSummary.coverageAllTruthCrossings=sum(coveredSyn)/max(1,synSummary.nTruthPlaneCrossings);
synSummary.classificationAccuracy=mean(strcmp({synSamples.truthStatus},{synSamples.fitStatus}));
synSummary.truthStatusNames = {'goal','ground','cross_out','timeout'};
synSummary.truthStatusCounts = [sum(strcmp({synSamples.truthStatus}', 'goal')), ...
    sum(strcmp({synSamples.truthStatus}', 'ground')), ...
    sum(strcmp({synSamples.truthStatus}', 'cross_out')), ...
    sum(strcmp({synSamples.truthStatus}', 'timeout'))];
synSummary.fitStatusNames = {'goal','ground','cross_out','timeout','fit_fail'};
synSummary.fitStatusCounts = [sum(strcmp({synSamples.fitStatus}', 'goal')), ...
    sum(strcmp({synSamples.fitStatus}', 'ground')), ...
    sum(strcmp({synSamples.fitStatus}', 'cross_out')), ...
    sum(strcmp({synSamples.fitStatus}', 'timeout')), ...
    sum(strcmp({synSamples.fitStatus}', 'fit_fail'))];

reach = isfinite([synSamples.yTrue]') & coverValidSyn;
spinNorm = sqrt(thTrue(:,5).^2 + thTrue(:,6).^2);
if any(reach)
    spinCut = median(spinNorm(reach));
else
    spinCut = NaN;
end
lowerMask = reach & spinNorm <= spinCut;
higherMask = reach & spinNorm > spinCut;
if any(lowerMask), lowerCov = sum(coveredSyn & lowerMask)/sum(lowerMask); else, lowerCov = NaN; end
if any(higherMask), higherCov = sum(coveredSyn & higherMask)/sum(higherMask); else, higherCov = NaN; end
synSummary.spinCoverage = struct('nReach', sum(reach), ...
    'medianSpin', spinCut, ...
    'lowerSpinDenom', sum(lowerMask), 'lowerSpinCoverage', lowerCov, ...
    'higherSpinDenom', sum(higherMask), 'higherSpinCoverage', higherCov);
synSummary.note = ['这是独立合成轨迹的数值覆盖检验，不是实验确认；覆盖分母只包含' ...
    '真值与反演都穿过球门平面且预测协方差有限的样本；另给包含失败的全分母覆盖率。'];

% ---------- 7. 模型偏差压力测试（每个情景重新拟合相同早期观测） ----------
scNames = {'bx=-1', 'bx=+1', 'wind_y=-3', 'wind_y=+3', ...
    'Cd exponent p=-0.5', 'Cd exponent p=+0.5', ...
    'b decay tau=0.5 s', 'b decay tau=1.0 s'};
scenarios = repmat(struct('name','', 'bx',0, 'wind',[0;0;0], ...
    'cdExponent',0, 'bDecayTau',inf, 'theta',nan(1,6), ...
    'rmse',NaN, 'residual',nan(1,6), 'status','', 'goal',false, ...
    'yGoal',NaN, 'zGoal',NaN, 'relativeDeviation',NaN, ...
    'boundaryHit',false, 'exitflag',0), numel(scNames), 1);
for k = 1:numel(scNames)
    cfgS = cfg;
    switch k
        case 1
            cfgS.bx = -1;
        case 2
            cfgS.bx = +1;
        case 3
            cfgS.w = [0; -3; 0];
        case 4
            cfgS.w = [0; +3; 0];
        case 5
            cfgS.cdExponent = -0.5;
        case 6
            cfgS.cdExponent = +0.5;
        case 7
            cfgS.bDecayTau = 0.5;
        case 8
            cfgS.bDecayTau = 1.0;
    end
    fitS = fitWithStarts(thetaHat, obsT, obs, cfgS);
    predS = predictGoal(fitS.theta, cfgS, 11.00);
    relDev = NaN;
    if predBase.goal && predS.goal
        relDev = norm([predS.yGoal; predS.zGoal] - ...
            [predBase.yGoal; predBase.zGoal]) / max(norm([predBase.yGoal; predBase.zGoal]), eps);
    end
    scenarios(k).name = scNames{k};
    scenarios(k).bx = cfgS.bx;
    scenarios(k).wind = cfgS.w;
    scenarios(k).cdExponent = cfgS.cdExponent;
    scenarios(k).bDecayTau = cfgS.bDecayTau;
    scenarios(k).theta = fitS.theta;
    scenarios(k).rmse = fitS.rmse;
    scenarios(k).residual = fitS.resid;
    scenarios(k).status = predS.status;
    scenarios(k).goal = predS.goal;
    scenarios(k).yGoal = predS.yGoal;
    scenarios(k).zGoal = predS.zGoal;
    scenarios(k).relativeDeviation = relDev;
    scenarios(k).boundaryHit = fitS.boundaryHit;
    scenarios(k).exitflag = fitS.exitflag;
    fprintf('[Mai_penalty] 模型情景 %d/%d：%s，状态=%s。\n', ...
        k, numel(scNames), scNames{k}, predS.status);
end

% ---------- 8. 正式图（英文轴/图例，SI单位，无标题） ----------
figFiles = {};
% 图1：轨迹与观测
fig = figure('Visible','off', 'Color','w');
hold on;
tGrid = linspace(0,predBase.tCross,120);
STraj = simulateAt(thetaHat, tGrid, cfg, 0.01);
plot3(STraj(:,1), STraj(:,2), STraj(:,3), 'b-', 'LineWidth', 2, ...
    'DisplayName','Model trajectory');
plot3(obs(:,1), obs(:,2), obs(:,3), 'ko', 'MarkerFaceColor','k', ...
    'MarkerSize',6, 'DisplayName','Observations');
plot3(11*ones(1,5), [-3.66 3.66 3.66 -3.66 -3.66], ...
    [0 0 2.44 2.44 0], 'r-', 'LineWidth', 1.5, 'DisplayName','Goal plane x = 11 m');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
grid on; view(3); legend('Location','best');
f1 = fullfile(imgDir, 'Mai罚球轨迹与观测.png');
exportgraphics(fig, f1, 'Resolution', 200); close(fig); figFiles{end+1} = f1;

% 图2：球门预测椭圆与Bootstrap散点
fig = figure('Visible','off', 'Color','w');
hold on;
idxBS = find(okBS & isfinite([bsSamples.yGoal]') & isfinite([bsSamples.zGoal]'));
if ~isempty(idxBS)
    scatter([bsSamples(idxBS).yGoal], [bsSamples(idxBS).zGoal], 12, ...
        [0.4 0.4 0.4], 'filled', 'DisplayName','Bootstrap predictions');
end
if predBase.goal && all(isfinite(C(:)))
    tt = linspace(0, 2*pi, 160);
    E = ell.vectors * diag(ell.axes) * [cos(tt); sin(tt)] + ...
        [predBase.yGoal; predBase.zGoal];
    plot(E(1,:), E(2,:), 'r-', 'LineWidth', 2, 'DisplayName','95% joint ellipse');
    plot(predBase.yGoal, predBase.zGoal, 'r+', 'MarkerSize',10, ...
        'LineWidth', 2, 'DisplayName','Plug-in prediction');
end
plot([-3.66 3.66 3.66 -3.66 -3.66], [0 0 2.44 2.44 0], 'k--', ...
    'DisplayName','Goal frame');
xlabel('y (m)'); ylabel('z (m)'); axis equal; grid on; legend('Location','best');
f2 = fullfile(imgDir, 'Mai球门预测椭圆与自助散点.png');
exportgraphics(fig, f2, 'Resolution', 200); close(fig); figFiles{end+1} = f2;

% 图3：噪声水平与边际区间半宽
fig = figure('Visible','off', 'Color','w');
HW = horzcat(noiseScenarios.marginal95HalfWidth_yz)';
plot(sigmas, HW(:,1), 'o-', 'LineWidth', 1.8, 'DisplayName','y half-width');
hold on;
plot(sigmas, HW(:,2), 's-', 'LineWidth', 1.8, 'DisplayName','z half-width');
xlabel('Assumed position noise \sigma (m)');
ylabel('95% marginal half-width (m)');
grid on; legend('Location','best');
f3 = fullfile(imgDir, 'Mai噪声水平与边际区间半宽.png');
exportgraphics(fig, f3, 'Resolution', 200); close(fig); figFiles{end+1} = f3;

% 图4：各单因素情景相对于基准的偏移，厘米单位便于辨认
fig = figure('Visible','off', 'Color','w');
hold on;
plot(0,0,'kp','MarkerFaceColor','k','MarkerSize',10,'DisplayName','Baseline');
for k=1:numel(scenarios)
    if all(isfinite([scenarios(k).yGoal,scenarios(k).zGoal]))
        plot(100*(scenarios(k).yGoal-predBase.yGoal),100*(scenarios(k).zGoal-predBase.zGoal), ...
            'o','MarkerSize',7,'DisplayName',scenarios(k).name);
    end
end
xlabel('Change in y (cm)'); ylabel('Change in z (cm)'); grid on; axis equal;
legend('Location','bestoutside','Interpreter','none');
f4=fullfile(imgDir,'Mai模型情景预测比较.png');
exportgraphics(fig,f4,'Resolution',200); close(fig); figFiles{end+1}=f4;

% ---------- 9. 汇总结果结构 ----------
result = struct();
result.mode = mode;
result.seed = seed;
result.software = struct('matlabVersion', version, 'computer', computer, ...
    'expectedRelease','R2024a', ...
    'toolboxes','Optimization Toolbox (lsqnonlin); Statistics and Machine Learning Toolbox (chi2inv, prctile)');
result.numericSettings = struct( ...
    'rk4Step_s', cfg.h, ...
    'rk4UsedForObservationFit', true, ...
    'predictionIntegrator','ode45', ...
    'predictionMaxStep_s', cfg.odeMaxStep, ...
    'predictionRelTol', 1e-9, 'predictionAbsTol', 1e-11, ...
    'maxFlightTime_s', 2.0, ...
    'maxTimeForStartState_s', 0.10, ...
    'lsqnonlinAlgorithm','trust-region-reflective', ...
    'lsqnonlinFunctionTolerance', 1e-12, ...
    'lsqnonlinStepTolerance', 1e-12, ...
    'lsqnonlinOptimalityTolerance', 1e-12, ...
    'lsqnonlinMaxFunctionEvaluations', 10000, ...
    'finiteDifferenceRelativeStep', 1e-6, ...
    'weightedResidual','(predicted - observed)/sigma', ...
    'sigma_m', cfg.sigma, ...
    'goalPlane_x_m', 11.00, ...
    'centerPassagePlane_x_m', 11.11, ...
    'goalGeometry', 'abs(y) < 3.66 - R, z < 2.44 - R, z >= R');
result.constants = struct('rho',cfg.rho, 'R',cfg.R, 'm',cfg.m, ...
    'K',cfg.K, 'g',cfg.g, 'r0',cfg.r0, ...
    'goalHalfWidth',3.66, 'goalHeight',2.44, 'sigmaAssumed',cfg.sigma);
result.observations = struct('times',obsT, 'positions',obs);
result.fit = fitBest;
result.fit.allStarts = fitStarts;
result.prediction = predBase;
result.prediction.centerPassage_x11_11 = predCenter;
result.uncertainty = unc;
result.noiseScenarios = noiseScenarios;
result.baselines = baselines;
result.bootstrap = struct('summary',bsSummary, 'samples',bsSamples);
result.verification = ver;
result.modelScenarios = scenarios;
result.synthetic = struct('summary',synSummary, 'samples',synSamples);
result.figures = figFiles;
result.timing = struct('elapsedSeconds', toc(tStart));
result.notes = { ...
    '本结果由MATLAB实际运行生成；仍未完成实验确认。'; ...
    ['bx固定为0，b=[0;by;bz]是有效Magnus系数向量，不能据此唯一恢复真实自旋。']; ...
    ['sigma=0.01 m 是独立各坐标高斯噪声的建模假设，不能由两次观测的6个数估计；' ...
     '另外给出0.002、0.005、0.02 m情景。']; ...
    ['r0=[0;0;0.11] 在拟合中暂按精确已知处理，这是明确限制。']; ...
    ['预测点是约束最大似然参数的plug-in点，不称为无条件唯一概率众数。']; ...
    ['联合95%椭圆和边际区间仅来自局部高斯线性化，不处理参数边界。']; ...
    ['Bootstrap的越过球门比例是重拟合稳定性频率，不是实际进球概率；' ...
     '失败与边界命中单独记录，未把筛选后的进球样本区间冒充无条件区间。']; ...
    ['模型偏差情景是逐个改变单一因素的压力测试，不是概率区间。']; ...
    ['合成反演覆盖检验不是实验确认。']; ...
    ['所有正式派生数据写入 data/派生数据，正式图片写入 img/。']};

% ---------- 10. 写出JSON/CSV/PNG/MAT ----------
jsonPath = fullfile(dataDir, 'Mai罚球结果.json');
try
    txt = jsonencode(result, 'PrettyPrint', true);
    fid = fopen(jsonPath, 'w', 'n', 'UTF-8');
    if fid < 0, error('无法打开JSON文件。'); end
    fprintf(fid, '%s', txt);
    fclose(fid);
catch ME
    warning('Mai_penalty:JsonWrite', 'JSON写出失败：%s', ME.message);
end

matPath = fullfile(dataDir, 'Mai罚球结果.mat');
save(matPath, 'result');

% 拟合残差CSV
Pfit = predictPositions(thetaHat, obsT, cfg);
Tfit = table(obsT(:), obs(:,1), obs(:,2), obs(:,3), ...
    Pfit(:,1), Pfit(:,2), Pfit(:,3), ...
    Pfit(:,1)-obs(:,1), Pfit(:,2)-obs(:,2), Pfit(:,3)-obs(:,3), ...
    'VariableNames', {'t_s','x_obs_m','y_obs_m','z_obs_m', ...
    'x_fit_m','y_fit_m','z_fit_m','rx_m','ry_m','rz_m'});
writetable(Tfit, fullfile(dataDir, 'Mai拟合残差.csv'));

% Bootstrap样本CSV
TH = vertcat(bsSamples.theta);
Tbs = table((1:nb)', TH(:,1), TH(:,2), TH(:,3), TH(:,4), TH(:,5), TH(:,6), ...
    [bsSamples.ok]', [bsSamples.exitflag]', [bsSamples.boundaryHit]', ...
    {bsSamples.status}', [bsSamples.goal]', [bsSamples.yGoal]', ...
    [bsSamples.zGoal]', [bsSamples.sse]', ...
    'VariableNames', {'index','vx','vy','vz','Cd','by','bz','ok','exitflag', ...
    'boundaryHit','status','goal','yGoal_m','zGoal_m','sse_m2'});
writetable(Tbs, fullfile(dataDir, 'Mai自助样本.csv'));

% 合成反演样本CSV
TT = vertcat(synSamples.thetaTrue);
THS = vertcat(synSamples.thetaHat);
Tsyn = table((1:nSyn)', TT(:,1), TT(:,2), TT(:,3), TT(:,4), TT(:,5), TT(:,6), ...
    THS(:,1), THS(:,2), THS(:,3), THS(:,4), THS(:,5), THS(:,6), ...
    [synSamples.ok]', [synSamples.boundaryHit]', [synSamples.truthGoal]', ...
    [synSamples.fitGoal]', [synSamples.yTrue]', [synSamples.zTrue]', ...
    [synSamples.yHat]', [synSamples.zHat]', [synSamples.coverValid]', ...
    [synSamples.covered]', {synSamples.truthStatus}', {synSamples.fitStatus}', ...
    'VariableNames', {'index','vx_true','vy_true','vz_true','Cd_true', ...
    'by_true','bz_true','vx_hat','vy_hat','vz_hat','Cd_hat','by_hat','bz_hat', ...
    'ok','boundaryHit','truthGoal','fitGoal','yTrue_m','zTrue_m', ...
    'yHat_m','zHat_m','coverValid','covered','truthStatus','fitStatus'});
writetable(Tsyn, fullfile(dataDir, 'Mai合成反演样本.csv'));

% 模型偏差情景CSV
Tsc = table({scenarios.name}', [scenarios.bx]', [scenarios.rmse]', ...
    [scenarios.exitflag]', [scenarios.boundaryHit]', {scenarios.status}', ...
    [scenarios.goal]', [scenarios.yGoal]', [scenarios.zGoal]', ...
    [scenarios.relativeDeviation]', ...
    'VariableNames', {'scenario','bx','rmse_m','exitflag','boundaryHit', ...
    'status','goal','yGoal_m','zGoal_m','relativeDeviation'});
writetable(Tsc, fullfile(dataDir, 'Mai模型偏差情景.csv'));

% 噪声情景CSV
Tns = table([noiseScenarios.sigma]', ...
    HW(:,1), HW(:,2), ...
    'VariableNames', {'sigma_m','yHalfWidth_m','zHalfWidth_m'});
writetable(Tns, fullfile(dataDir, 'Mai噪声情景.csv'));

fprintf('[Mai_penalty] 完成。JSON=%s，耗时=%.2f s。\n', jsonPath, result.timing.elapsedSeconds);
end

% =====================================================================
% 局部函数
% =====================================================================
function [best, allStarts] = fitWithStarts(starts, times, obs, cfg)
%FITWITHSTARTS 有界加权最小二乘多初值拟合
nS = size(starts, 1);
opts = optimoptions('lsqnonlin', 'Display','off', ...
    'Algorithm','trust-region-reflective', ...
    'FunctionTolerance',1e-12, 'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-12, ...
    'MaxFunctionEvaluations',10000, 'MaxIterations',2000);
allStarts = repmat(struct('theta',nan(1,6), 'sse',NaN, 'rmse',NaN, ...
    'exitflag',0, 'boundaryHit',false), nS, 1);
best = struct('theta',nan(1,6), 'sse',inf, 'rmse',NaN, ...
    'resid',nan(1,6), 'exitflag',0, 'boundaryHit',false);
for k = 1:nS
    % 初值夹到上下界内部微小余量，避免TRR算法从边界上启动
    x0 = min(max(starts(k,:), cfg.lb + 1e-9), cfg.ub - 1e-9);
    x0f = x0(cfg.freeMask);
    [xf, ~, ~, exitflag] = lsqnonlin( ...
        @(x) fitResidual(x, times, obs, cfg), ...
        x0f, cfg.lb(cfg.freeMask), cfg.ub(cfg.freeMask), opts);
    theta = cfg.thetaFixed;
    theta(cfg.freeMask) = xf;
    P = predictPositions(theta, times, cfg);
    resid = P - obs;
    sse = sum(resid(:).^2);
    rmse = sqrt(mean(resid(:).^2));
    hit = hitBoundary(theta, cfg);
    allStarts(k).theta = theta;
    allStarts(k).sse = sse;
    allStarts(k).rmse = rmse;
    allStarts(k).exitflag = exitflag;
    allStarts(k).boundaryHit = hit;
    if exitflag>0 && sse < best.sse
        best.theta = theta;
        best.sse = sse;
        best.rmse = rmse;
        best.resid = reshape(resid, 1, []);
        best.exitflag = exitflag;
        best.boundaryHit = hit;
    end
end
end

function r = fitResidual(x, times, obs, cfg)
theta = cfg.thetaFixed;
theta(cfg.freeMask) = x;
P = predictPositions(theta, times, cfg);
r = (P - obs) ./ cfg.sigma;
r = r(:);
end

function tf = hitBoundary(theta, cfg)
tf = false;
for j = find(cfg.freeMask)
    tol = 1e-4 * max(1, abs(cfg.lb(j)));
    if abs(theta(j) - cfg.lb(j)) <= tol || abs(theta(j) - cfg.ub(j)) <= tol
        tf = true;
    end
end
end

function P = predictPositions(theta, times, cfg)
[P, ~] = rk4State(theta, times, cfg);
end

function [P, stateEnd] = rk4State(theta, times, cfg)
%RK4STATE 固定步长RK4，精确记录到给定观测时刻
times = times(:)';
assert(all(diff(times)>=0) && all(times>=0),'Invalid observation times');
state = [cfg.r0; theta(1:3)'];
P = nan(numel(times),3);
t = 0;
for j = 1:numel(times)
    while t < times(j)-1e-13
        h = min(cfg.h,times(j)-t);
        k1=stateDeriv(t,state,cfg,theta);
        k2=stateDeriv(t+h/2,state+h*k1/2,cfg,theta);
        k3=stateDeriv(t+h/2,state+h*k2/2,cfg,theta);
        k4=stateDeriv(t+h,state+h*k3,cfg,theta);
        state=state+h*(k1+2*k2+2*k3+k4)/6;
        t=t+h;
    end
    P(j,:)=state(1:3)';
end
stateEnd=state;
end

function ds = stateDeriv(t, s, cfg, theta)
r = s(1:3); %#ok<NASGU>
v = s(4:6);
u = v - cfg.w;
Cd = theta(4);
if cfg.cdExponent ~= 0
    Cd = Cd * (norm(u)/cfg.cdRefSpeed)^cfg.cdExponent;
end
b = [cfg.bx; theta(5); theta(6)];
if isfinite(cfg.bDecayTau) && cfg.bDecayTau > 0
    b = b * exp(-t/cfg.bDecayTau);
end
a = cfg.g - cfg.K*Cd*norm(u)*u + cross(b, u);
ds = [v; a];
end

function S = simulateAt(theta, times, cfg, maxStep)
opts = odeset('RelTol',1e-11, 'AbsTol',1e-13, 'MaxStep',maxStep);
sol = ode45(@(t,s) stateDeriv(t,s,cfg,theta), ...
    [0,max(times)], [cfg.r0;theta(1:3)'],opts);
S = deval(sol,times)';
end

function pred = predictGoal(theta, cfg, targetX)
%PREDICTGOAL 从t=0.1状态开始事件积分，定位x=targetX或触地z=R
[~, s01] = rk4State(theta, [0.05, 0.10], cfg);
opts = odeset('RelTol',1e-9, 'AbsTol',1e-11, ...
    'MaxStep',cfg.odeMaxStep, ...
    'Events', @(t,s) goalEvents(t, s, cfg, targetX));
[~, Y, te, ye, ie] = ode45(@(t,s) stateDeriv(t,s,cfg,theta), ...
    [0.10, 2.0], s01, opts); %#ok<ASGLU>
pred = struct('status','timeout', 'goal',false, ...
    'yGoal',NaN, 'zGoal',NaN, 'tCross',NaN, 'point',[NaN;NaN], ...
    'targetX',targetX);
if isempty(te) || isempty(ie)
    return;
end
first = ie(1);
ye1 = ye(1,:);
pred.tCross = te(1);
if first == 2
    pred.status = 'ground';
    return;
end
y = ye1(2);
z = ye1(3);
pred.yGoal=y; pred.zGoal=z; pred.point=[y;z];
if abs(y) < 3.66 - cfg.R && z < 2.44 - cfg.R && z >= cfg.R
    pred.status = 'goal';
    pred.goal = true;
    pred.yGoal = y;
    pred.zGoal = z;
    pred.point = [y; z];
else
    pred.status = 'cross_out';
end
end

function [value, isterminal, direction] = goalEvents(~, s, cfg, targetX)
value = [s(1) - targetX; s(3) - cfg.R];
isterminal = [1; 1];
direction = [1; -1];
end

function [J, info] = obsJacobian(theta, cfg, times)
base = predictPositions(theta, times, cfg);
base = base(:);
J = nan(numel(base), 6);
dAll = zeros(1,6);
for j = 1:6
    d = 1e-6 * max(1, abs(theta(j)));
    dAll(j) = d;
    tp = theta; tm = theta;
    tp(j) = tp(j) + d;
    tm(j) = tm(j) - d;
    pValid = tp(j) <= cfg.ub(j) + eps;
    mValid = tm(j) >= cfg.lb(j) - eps;
    if pValid, pp = predictPositions(tp, times, cfg); pp = pp(:); end
    if mValid, pm = predictPositions(tm, times, cfg); pm = pm(:); end
    if pValid && mValid
        J(:,j) = (pp - pm) / (tp(j) - tm(j));
    elseif pValid
        J(:,j) = (pp - base) / d;
    elseif mValid
        J(:,j) = (base - pm) / d;
    end
end
info = struct('perturbation', dAll);
end

function [G, info] = goalJacobian(theta, cfg)
pb = predictGoal(theta, cfg, 11.00);
G = nan(2, 6);
dAll = zeros(1,6);
if ~all(isfinite([pb.yGoal,pb.zGoal]))
    info = struct('perturbation', dAll);
    return;
end
gb = [pb.yGoal; pb.zGoal];
for j = 1:6
    d = 1e-6 * max(1, abs(theta(j)));
    dAll(j) = d;
    tp = theta; tm = theta;
    tp(j) = tp(j) + d;
    tm(j) = tm(j) - d;
    pValid = tp(j) <= cfg.ub(j) + eps;
    mValid = tm(j) >= cfg.lb(j) - eps;
    if pValid
        pp = predictGoal(tp, cfg, 11.00);
    else
        pp = struct('goal',false, 'yGoal',NaN, 'zGoal',NaN);
    end
    if mValid
        pm = predictGoal(tm, cfg, 11.00);
    else
        pm = struct('goal',false, 'yGoal',NaN, 'zGoal',NaN);
    end
    if pValid && mValid && all(isfinite([pp.yGoal,pp.zGoal,pm.yGoal,pm.zGoal]))
        G(:,j) = ([pp.yGoal; pp.zGoal] - [pm.yGoal; pm.zGoal]) / (tp(j) - tm(j));
    elseif pValid && all(isfinite([pp.yGoal,pp.zGoal]))
        G(:,j) = ([pp.yGoal; pp.zGoal] - gb) / d;
    elseif mValid && all(isfinite([pm.yGoal,pm.zGoal]))
        G(:,j) = (gb - [pm.yGoal; pm.zGoal]) / d;
    end
end
info = struct('perturbation', dAll);
end

function e = ellipseFromCov(C)
e = struct('axes',[NaN;NaN], 'vectors',nan(2,2));
if any(~isfinite(C(:)))
    return;
end
C = (C + C')/2;
[V, D] = eig(C);
[dd, ord] = sort(diag(D), 'descend');
V = V(:, ord);
e.axes = sqrt(chi2inv(0.95, 2) * max(dd, 0));
e.vectors = V;
end

function p = analyticPrediction(r0, v0, a, cfg)
%ANALYTICPREDICTION 常加速度解析预测x=11最早正实根
tx = positiveQuadraticRoots(0.5*a(1), v0(1), -11);
tg = positiveQuadraticRoots(0.5*a(3), v0(3), r0(3) - cfg.R);
p = struct('status','no_goal', 'goal',false, 'yGoal',NaN, ...
    'zGoal',NaN, 'tCross',NaN, 'point',[NaN;NaN]);
if isempty(tx)
    return;
end
if ~isempty(tg) && tg(1) < tx(1)
    p.status = 'ground';
    p.tCross = tg(1);
    return;
end
t = tx(1);
pos = r0 + v0*t + 0.5*a*t^2;
p.tCross = t;
p.yGoal=pos(2); p.zGoal=pos(3); p.point=pos(2:3)';
if abs(pos(2)) < 3.66 - cfg.R && pos(3) < 2.44 - cfg.R && pos(3) >= cfg.R
    p.status = 'goal';
    p.goal = true;
    p.yGoal = pos(2);
    p.zGoal = pos(3);
    p.point = [pos(2); pos(3)];
else
    p.status = 'cross_out';
end
end

function ts = positiveQuadraticRoots(a, b, c)
if abs(a) < 1e-14
    if abs(b) < 1e-14
        ts = [];
    else
        ts = -c/b;
        ts = sort(ts(ts > 0));
    end
    return;
end
disc = b^2 - 4*a*c;
if disc < 0
    ts = [];
else
    rr = [(-b + sqrt(disc))/(2*a), (-b - sqrt(disc))/(2*a)];
    ts = sort(rr(rr > 0));
end
end

function b = makeBaseline(name, theta, v0, accel, P, obs, pred)
resid = P - obs;
b = struct('name',name, 'theta',theta, 'v0',v0, 'accel',accel, ...
    'rmse',sqrt(mean(resid(:).^2)), 'residual',reshape(resid,1,[]), ...
    'status',pred.status, 'goal',pred.goal, ...
    'yGoal',pred.yGoal, 'zGoal',pred.zGoal, 'tCross',pred.tCross);
end
