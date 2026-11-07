function out=penalty_scenario(action,varargin)
%PENALTY_SCENARIO 全程参数约束反演；不再拼接t=.1前后的不同物理。
switch lower(action)
 case 'defaults'
  out=penalty_inverse('defaults');d=loadData();out.Cd=d.fit.theta(4);
 case 'start'
  p=varargin{1};d=loadData();f=penalty_inverse('fit',p,d.observations.positions);
  sim=penalty_inverse('propagate',f);sim.base.observations=d.observations;
  refp=penalty_inverse('defaults');refp.Cd=d.fit.theta(4);
  sim.reference=reference(refp,d.observations.positions);out=sim;
 case 'sample'
  sim=varargin{1};t=varargin{2};
  assert(isscalar(t)&&isfinite(t)&&t>=0&&t<=sim.tEnd+1e-9,'时间越界。');
  s=deval(sim.sol,min(t,sim.tEnd));[~,force]=ball_rhs(t,s,sim.theta,sim.cfg);
  out=struct('state',s,'forces',force,'wind',force.wind,'Cd',force.Cd);
 case 'curve'
  sim=varargin{1};t=unique([linspace(0,sim.tEnd,240) .05 .1]);t=t(t<=sim.tEnd);
  out=struct('t',t(:),'state',deval(sim.sol,t)');
 case 'family'
  d=loadData();out=struct([]);
  if isfield(d,'trajectoryFamily')
   a=d.trajectoryFamily.members;out=a([a.maxResidual]<=varargin{1}&[a.exitflag]>0&[a.physicalAdmissible]);
   for k=1:numel(out),out(k).accepted=true;end
  end
 case 'branch'
  error('penalty_scenario:RefitRequired','环境参数必须从t=0重新反演；使用start操作。');
 otherwise,error('未知仿真操作：%s',action);
end
end

function d=loadData()
persistent cache stamp
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
p=fullfile(root,'data','派生数据','罚球结果.json');m=dir(p);
assert(~isempty(m),'请先运行 penalty(''full'') 生成正式结果。');tag=[m.datenum m.bytes];
if isempty(cache)||~isequal(stamp,tag),cache=jsondecode(fileread(p));stamp=tag;end
d=cache;
end

function ref=reference(p,obs)
persistent cached signature
key=[p.Cd obs(:)'];
if isempty(cached)||~isequal(key,signature)
 sim=penalty_inverse('propagate',penalty_inverse('fit',p,obs));t=linspace(0,sim.tEnd,240);
 cached=struct('t',t(:),'state',deval(sim.sol,t)','endpoint',sim.endpoint,'tEnd',sim.tEnd);signature=key;
end
ref=cached;
end
