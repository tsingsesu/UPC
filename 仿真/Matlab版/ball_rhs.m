function [ds, forces] = ball_rhs(t, s, theta, cfg)
%BALL_RHS 点球飞行右端项（与 penalty.stateDeriv 逐字一致）
%
%   [ds, forces] = ball_rhs(t, s, theta, cfg)
%
%   状态 s = [x;y;z; vx;vy;vz]，theta = [vx vy vz Cd by bz]。
%   ds = [v; g - K*Cd*||u||*u + cross(b,u)]，u = v - w。
%   forces 给出重力/阻力/Magnus 三维力向量（N），质量取 cfg.m：
%       Fg = m*g, Fd = -m*K*Cd*||u||*u, Fm = m*cross(b,u)。
%   cfg 支持 g、K、w、bx、cdExponent、cdRefSpeed、bDecayTau。
%   b = [bx; by; bz] 是有效Magnus系数向量，单位 1/s。

s = s(:);
v = s(4:6);
wind = cfg.w(:); effectiveK = cfg.K;
CdMultiplier = 1; spinMultiplier = 1; buoyancy = zeros(3,1);
if isfield(cfg,'environment') && t >= cfg.activationTime
    e = cfg.environment;
    wind = [e.windX;e.windY;e.windZ];
    phase = (t-e.gustStart)/e.gustDuration;
    if phase>0 && phase<1
        wind = wind + sin(pi*phase)^2*[0;e.gustY;e.gustZ];
    end
    % 确定性平滑扰动，用于压力测试，不是实测湍流概率模型。
    age=max(0,t-cfg.physicsEpoch); ramp=min(1,age/.02);
    wind=wind+e.turbulence*ramp*[sin(2*pi*7*age);sin(2*pi*11*age+.7);.5*sin(2*pi*5*age+1.3)];
    effectiveK=e.density*pi*cfg.R^2/(2*cfg.m);
    CdMultiplier=e.CdScale;
    spinMultiplier=e.spinScale;
    if isfield(cfg,'referenceMass')
        spinMultiplier=spinMultiplier*(e.density/cfg.referenceDensity) ...
            *(cfg.R/cfg.referenceRadius)^3*(cfg.referenceMass/cfg.m);
    end
    if e.spinTau>0, spinMultiplier=spinMultiplier*exp(-max(0,t-cfg.physicsEpoch)/e.spinTau); end
    % 空气浮力向上；默认禁用，确保旧基线数值完全兼容。
    if isfield(e,'buoyancy') && e.buoyancy
        buoyancy=-e.density*(4*pi*cfg.R^3/3)*cfg.g(:);
    end
end
u = v - wind;
Cd = theta(4)*CdMultiplier;
exponent=cfg.cdExponent;
if isfield(cfg,'environment') && t>=cfg.activationTime, exponent=cfg.environment.cdExponent; end
if exponent ~= 0
    Cd=Cd*(max(norm(u),1e-9)/cfg.cdRefSpeed)^exponent;
end
b = [cfg.bx; theta(5); theta(6)];
if isfield(cfg,'environment') && t>=cfg.activationTime
    b(1)=cfg.environment.bx;
end
b=b*spinMultiplier;
if isfinite(cfg.bDecayTau) && cfg.bDecayTau > 0
    b = b * exp(-t/cfg.bDecayTau);
end
a = cfg.g(:) - effectiveK*Cd*norm(u)*u + cross(b, u) + buoyancy/cfg.m;
ds = [v; a];

% 与加速度表达式严格对应的力，仅按质量换算为牛顿
forces = struct( ...
    'gravity', cfg.m*cfg.g(:), ...
    'drag',    -cfg.m*effectiveK*Cd*norm(u)*u, ...
    'magnus',  cfg.m*cross(b, u), 'buoyancy',buoyancy, ...
    'wind',wind,'Cd',Cd,'relativeVelocity',u);
end
