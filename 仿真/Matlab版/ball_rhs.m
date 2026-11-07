function [ds, forces] = ball_rhs(t, s, theta, cfg)
%BALL_RHS 全程生效的重力、相对气流平方阻力和有效Magnus力。
% 所有参数从出脚t=0起使用，无0.10秒激活门槛。
% theta=[v0x v0y v0z Cd by bz]，b单位s^-1，不等同实际角速度。
% 可选Cd速度指数及b衰减供既有单因素压力测试使用。
s=s(:);v=s(4:6);u=v-cfg.w(:);Cd=theta(4);
if cfg.cdExponent~=0
 Cd=Cd*(max(norm(u),1e-9)/cfg.cdRefSpeed)^cfg.cdExponent;
end
b=[cfg.bx;theta(5);theta(6)];
if isfinite(cfg.bDecayTau)&&cfg.bDecayTau>0,b=b*exp(-t/cfg.bDecayTau);end
Fd=-cfg.m*cfg.K*Cd*norm(u)*u;Fm=cfg.m*cross(b,u);Fg=cfg.m*cfg.g(:);
ds=[v;(Fg+Fd+Fm)/cfg.m];
if nargout>1
 forces=struct('gravity',Fg,'drag',Fd,'magnus',Fm,'buoyancy',zeros(3,1), ...
  'wind',cfg.w(:),'Cd',Cd,'relativeVelocity',u);
end
end
