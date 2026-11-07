function out = penalty_scenario(action,varargin)
%PENALTY_SCENARIO 因果分段轨迹：已飞行的状态不因未来调参被改写。
switch lower(action)
    case 'defaults'
        out=struct('windX',0,'windY',0,'windZ',0,'gustY',0,'gustZ',0, ...
            'gustStart',.16,'gustDuration',.16,'turbulence',0, ...
            'CdScale',1,'cdExponent',0,'spinScale',1,'bx',0,'spinTau',0, ...
            'density',1.225,'mode','conditioned','speedScale',1, ...
            'yawDeg',0,'elevationDeg',0,'mass',.43,'radius',.11,'buoyancy',false,'maxStep',.002);
    case 'start'
        p=validateParams(varargin{1}); base=loadBase();
        sim=struct('base',base,'reference',base.reference,'params',p, ...
            'segments',[],'edits',struct('time',{},'params',{}));
        cfg=base.cfg; theta=base.theta;
        if strcmp(p.mode,'conditioned')
            % 球属性和起脚条件在此模式由基准约束，不在飞行中更换足球。
            p.mass=cfg.m; p.radius=cfg.R; p.speedScale=1;p.yawDeg=0;p.elevationDeg=0;
            sim.params=p;
            seg=makeSegment(0,.1,[cfg.r0;theta(1:3)'],theta,cfg,false);
            sim.segments=seg;
            s=deval(seg.sol,.1); t0=.1;
        else
            cfg.m=p.mass;cfg.R=p.radius;cfg.r0(3)=p.radius;
            cfg.K=cfg.rho*pi*cfg.R^2/(2*cfg.m);
            v=theta(1:3);speed=norm(v)*p.speedScale;
            yaw=atan2(v(2),v(1))+deg2rad(p.yawDeg);
            pitch=atan2(v(3),norm(v(1:2)))+deg2rad(p.elevationDeg);
            theta(1:3)=speed*[cos(pitch)*cos(yaw),cos(pitch)*sin(yaw),sin(pitch)];
            s=[cfg.r0;theta(1:3)'];t0=0;
        end
        sim.initialSpeed=norm(theta(1:3));
        cfg.environment=p;cfg.activationTime=t0;cfg.physicsEpoch=t0;
        future=makeSegment(t0,2,s,theta,cfg,true);
        sim.segments=[sim.segments future];
        sim.edits=struct('time',t0,'params',p);
        out=finish(sim,future);
    case 'branch'
        sim=varargin{1};t=varargin{2};p=validateParams(varargin{3});
        immutable={'mode','speedScale','yawDeg','elevationDeg','mass','radius'};
        for k=1:numel(immutable)
            f=immutable{k};assert(isequal(p.(f),sim.params.(f)), ...
                'penalty_scenario:RestartRequired','起脚条件和球属性只能在重开球时调整。');
        end
        assert(isfinite(t)&&t>=0&&t<=sim.tEnd+1e-9,'无效的分支时刻。');
        if strcmp(p.mode,'conditioned'), t=max(.1,t); end
        assert(t<sim.tEnd-1e-9,'已经终止，请重播或拖回较早时刻再调整。');
        now=sampleAt(sim,t);s=now.state;
        idx=find([sim.segments.t0]<=t+1e-12 & [sim.segments.tf]>=t-1e-12,1,'last');
        old=sim.segments(idx);cfg=old.cfg;theta=old.theta;
        keep=sim.segments([sim.segments.t0]<t-1e-12);
        if ~isempty(keep),keep(end).tf=t;end
        cfg.environment=p;cfg.activationTime=t;
        % 自旋衰减的时间原点在起飞/最后观测时刻，不随每次改风速重置。
        if strcmp(p.mode,'conditioned'),cfg.physicsEpoch=.1;else,cfg.physicsEpoch=0;end
        future=makeSegment(t,2,s,theta,cfg,true);
        sim.segments=[keep future];sim.params=p;
        if ~isempty(sim.edits),sim.edits=sim.edits([sim.edits.time]<t-1e-12);end
        sim.edits(end+1)=struct('time',t,'params',p);
        out=finish(sim,future);
    case 'sample'
        out=sampleAt(varargin{1},varargin{2});
    case 'curve'
        sim=varargin{1};tt=[];ss=[];
        for k=1:numel(sim.segments)
            q=sim.segments(k);ts=linspace(q.t0,q.tf,max(2,ceil((q.tf-q.t0)/.002)+1));
            tt=[tt;ts(:)];ss=[ss;deval(q.sol,ts)']; %#ok<AGROW>
        end
        [tt,ix]=unique(tt,'last');out=struct('t',tt,'state',ss(ix,:));
    otherwise
        error('未知的实时仿真操作：%s',action);
end
end

function base=loadBase()
% 缓存仅用于一个MATLAB进程；文件修改后通过时间/大小自动失效。
persistent cache stamp
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
path=fullfile(root,'data','派生数据','罚球结果.json');meta=dir(path);
assert(~isempty(meta),'缺少正式基准结果。');tag=[meta.datenum,meta.bytes];
if isempty(cache)||~isequal(tag,stamp)
    d=jsondecode(fileread(path));c=d.constants;
    cfg=struct('R',c.R,'m',c.m,'rho',c.rho,'g',c.g(:),'r0',c.r0(:), ...
        'K',c.K,'referenceMass',c.m,'referenceRadius',c.R,'referenceDensity',c.rho,'w',[0;0;0],'bx',0,'cdExponent',0,'cdRefSpeed',30,'bDecayTau',inf);
    theta=d.fit.theta(:)';
    refseg=makeSegment(0,2,[cfg.r0;theta(1:3)'],theta,cfg,true);
    ts=unique([linspace(0,refseg.tf,240),d.observations.times(:)']);
    ref=struct('t',ts(:),'state',deval(refseg.sol,ts)', ...
        'tEnd',refseg.tf,'endpoint',deval(refseg.sol,refseg.tf));
    assert(norm(ref.endpoint(2:3)-[d.prediction.yGoal;d.prediction.zGoal])<1e-7, ...
        '参考轨迹不匹配正式计算。');
    cache=struct('cfg',cfg,'theta',theta,'observations',d.observations,'reference',ref);
    stamp=tag;
end
base=cache;
end

function seg=makeSegment(t0,tf,s,theta,cfg,events)
assert(tf>t0&&all(isfinite(s)),'无效的积分初始条件。');
step=.002;
if isfield(cfg,'environment'),step=cfg.environment.maxStep;end
opts=odeset('RelTol',1e-9,'AbsTol',1e-11,'MaxStep',step);
if events
    opts=odeset(opts,'Events',@(t,y)terminal(t,y,t0,cfg));
end
sol=ode45(@(t,y)ball_rhs(t,y,theta,cfg),[t0 tf],s,opts);
assert(all(isfinite(sol.y(:))),'积分出现非有限状态。');
assert(sol.x(end)>t0,'积分未推进。');
seg=struct('t0',t0,'tf',sol.x(end),'sol',sol,'cfg',cfg,'theta',theta,'status','timeout');
if events&&~isempty(sol.ie)
    if sol.ie(end)==1,seg.status='plane';else,seg.status='ground';end
end
end

function [value,stop,direction]=terminal(t,s,t0,cfg)
value=[s(1)-11;s(3)-cfg.R];
if t<=t0+1e-9 && t0==0,value(2)=max(value(2),1e-8);end
stop=[1;1];direction=[1;-1];
end

function sim=finish(sim,last)
sim.tEnd=last.tf;sim.endpoint=deval(last.sol,last.tf);sim.status=last.status;
sim.isOpening=strcmp(sim.status,'plane')&&abs(sim.endpoint(2))<3.66-sim.params.radius ...
    &&sim.endpoint(3)>=sim.params.radius&&sim.endpoint(3)<2.44-sim.params.radius;
end

function out=sampleAt(sim,t)
assert(isscalar(t)&&isfinite(t)&&t>=-1e-10&&t<=sim.tEnd+1e-9,'查询时间超出轨迹。');
t=min(max(0,t),sim.tEnd);
idx=find([sim.segments.t0]<=t+1e-11 & [sim.segments.tf]>=t-1e-11,1,'last');
assert(~isempty(idx),'分段轨迹缺失。');q=sim.segments(idx);
s=deval(q.sol,min(max(t,q.t0),q.tf));[~,forces]=ball_rhs(t,s,q.theta,q.cfg);
params=sim.params;if isfield(q.cfg,'environment'),params=q.cfg.environment;end
out=struct('state',s,'forces',forces,'wind',forces.wind,'Cd',forces.Cd,'params',params);
end

function p=validateParams(p)
d=penalty_scenario('defaults');names=fieldnames(d);
for k=1:numel(names)
    f=names{k};if ~isfield(p,f),p.(f)=d.(f);end
    if strcmp(f,'mode'),continue;end
    assert(isscalar(p.(f))&&isreal(p.(f))&&isfinite(p.(f)),['参数无效：' f]);
end
assert(any(strcmp(p.mode,{'conditioned','free'})),'未知起脚模式。');
assert(p.gustDuration>0&&p.gustStart>=0&&p.CdScale>0&&p.density>0 ...
    &&p.maxStep>0&&p.mass>0&&p.radius>0&&p.spinTau>=0&&p.speedScale>0&&p.turbulence>=0, ...
    '参数必须在物理可行范围内。');
end
