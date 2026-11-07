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
u = v - cfg.w(:);
Cd = theta(4);
if cfg.cdExponent ~= 0
    Cd = Cd * (norm(u)/cfg.cdRefSpeed)^cfg.cdExponent;
end
b = [cfg.bx; theta(5); theta(6)];
if isfinite(cfg.bDecayTau) && cfg.bDecayTau > 0
    b = b * exp(-t/cfg.bDecayTau);
end
a = cfg.g(:) - cfg.K*Cd*norm(u)*u + cross(b, u);
ds = [v; a];

% 与加速度表达式严格对应的力，仅按质量换算为牛顿
forces = struct( ...
    'gravity', cfg.m*cfg.g(:), ...
    'drag',    -cfg.m*cfg.K*Cd*norm(u)*u, ...
    'magnus',  cfg.m*cross(b, u));
end
