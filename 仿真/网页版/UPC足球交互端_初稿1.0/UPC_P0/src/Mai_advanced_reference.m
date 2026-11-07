function reference=Mai_advanced_reference()
%MAI_ADVANCED_REFERENCE 独立ode45实现M4，并输出浏览器跨实现参考。
% 运行：addpath('src'); Mai_advanced_reference();
% 方程依据：https://arxiv.org/html/1710.02784，式5/7/9及表1。
% 文献非旋转Cd数据拟合；Teamgeist升力含物理插值。线性CL、斜轴推广、
% 球壳惯量和衰减力矩明确为情景假设。不修改M0—M3主程序或团队正式结果。
c=struct('schemaVersion','1.0','model','M4','mode','forward', ...
    'initial',[0,0,.11],'velocity',[30,5,6],'omega',[0,0,50], ...
    'quaternion',[1,0,0,0],'m',.43,'R',.11,'rho',1.225,'g',9.80665, ...
    'wind',[0,0,0],'windShear',[0,0,0],'gustAmplitude',[0,0,0],'gustPeriod',1, ...
    'sigma',.01,'h',.001,'fitH',.005,'horizon',2,'goalX',11, ...
    'goalHalfWidth',3.66,'goalHeight',2.44,'Cd',.3,'b',[0,0,0], ...
    'constantAcceleration',[-16,-16,-12],'spinDecayTau','infinite', ...
    'kappa',2/3,'dragModel','constant','ballType','Teamgeist', ...
    'liftModel','linear','liftSlope',.8,'mu',1.81e-5);
c.observations=[struct('t',.05,'position',[1.5,.45,.44]),struct('t',.1,'position',[2.96,.86,.74])];
names={'constant_linear','teamgeist2018','brazuca_nonspinning_drag_linear_lift', ...
    'decay_shear_gust','outside_initial_domain'};
configs=repmat(c,1,5);
configs(2).dragModel='paper2018';configs(2).liftModel='teamgeist2018';
configs(3).dragModel='paper2018';configs(3).ballType='Brazuca';
configs(4).spinDecayTau=.5;configs(4).wind=[1,2,0];
configs(4).windShear=[.2,.3,0];configs(4).gustAmplitude=[0,1,.5];configs(4).gustPeriod=.7;
configs(5).dragModel='paper2018';configs(5).velocity=[36,0,1];
rows=cell(1,numel(names));
for k=1:numel(names)
    cc=configs(k);s0=[cc.initial,cc.velocity,cc.omega,cc.quaternion]';
    p=struct('status','no_crossing_within_horizon','eventTime',cc.horizon, ...
        'tGoal',NaN,'yGoal',NaN,'zGoal',NaN);
    sampleTimes=[0,.05,.1,.2,.3];
    if domainMargin(0,s0,cc)<0
        p.status='outside_model_domain';p.eventTime=0;
        sampleTimes=0;states=s0';finalState=s0';
    else
        opts=odeset('RelTol',1e-11,'AbsTol',1e-12,'MaxStep',.002, ...
            'Events',@(t,s) events(t,s,cc));
        sol=ode45(@(t,s) rhs(t,s,cc),[0,cc.horizon],s0,opts);
        finalState=sol.y(:,end)';
        if isfield(sol,'xe')&&~isempty(sol.xe)
            eventIndex=sol.ie(1);p.eventTime=sol.xe(1);finalState=sol.ye(:,1)';
            if eventIndex==1
                p.status='crossed_goal_plane';p.tGoal=p.eventTime;
                p.yGoal=finalState(2);p.zGoal=finalState(3);
            elseif eventIndex==2,p.status='ground_before_goal';
            else,p.status='outside_model_domain';end
        end
        sampleTimes=sampleTimes(sampleTimes<=p.eventTime);
        states=deval(sol,sampleTimes)';
    end
    assert(all(isfinite(states),'all'));
    rows{k}=struct('id',names{k},'config',cc,'prediction',p, ...
        'sampleTimes',sampleTimes,'states',states,'finalState',finalState, ...
        'quaternionMaxNormError',max(abs(sqrt(sum(states(:,10:13).^2,2))-1)));
end
reference=struct('schemaVersion','1.0','matlabVersion',version, ...
    'method','independently encoded 13-state ode45; rtol 1e-11, atol 1e-12, maxStep .002s', ...
    'stateOrder',{{'x','y','z','vx','vy','vz','omega_x','omega_y','omega_z','q_w','q_x','q_y','q_z'}}, ...
    'source','https://arxiv.org/html/1710.02784','cases',[rows{:}], ...
    'scope','跨实现数值一致性；不是现实实验确认；线性升力和力矩为假设。', ...
    'tolerances',struct('position_m',1e-4,'time_s',1e-5,'omega_rad_s',1e-5,'quaternion',1e-5));
targetDir=fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','验证基准');
if ~isfolder(targetDir),mkdir(targetDir);end
fid=fopen(fullfile(targetDir,'M4跨实现对照.json'),'w','n','UTF-8');
if fid<0,error('无法写出M4参考。');end
w=fprintf(fid,'%s',jsonencode(reference,'PrettyPrint',true));closed=fclose(fid);
assert(w>0&&closed==0,'M4参考写出失败。');
fprintf('M4独立MATLAB ode45：%d个对照算例已生成。\n',numel(rows));
end

function ds=rhs(t,s,c)
[drag,lift]=forces(t,s,c);
omega=s(7:9);q=s(10:13);
if ischar(c.spinDecayTau),omegaDot=zeros(3,1);else,omegaDot=-omega/c.spinDecayTau;end
% world-frame angular velocity, scalar-first body-to-world quaternion.
qdot=.5*[-dot(omega,q(2:4));q(1)*omega+cross(omega,q(2:4))];
ds=[s(4:6);[0;0;-c.g]+(drag+lift)/c.m;omegaDot;qdot];
end

function [drag,lift,U,S]=forces(t,s,c)
air=c.wind'+c.windShear'*s(3)+c.gustAmplitude'*sin(2*pi*t/c.gustPeriod);
w=s(4:6)-air;U=norm(w);omega=s(7:9);om=norm(omega);
drag=zeros(3,1);lift=zeros(3,1);S=NaN;
if U==0,return;end
S=c.R*om/U;d=c.Cd;
if strcmp(c.dragModel,'paper2018')
    switch c.ballType
        case 'Tango12',a=.5452;vc=12.86;vs=1.304;bmin=.1657;bmax=.1953;vmin=16.22;
        case 'Teamgeist',a=.4927;vc=12.58;vs=1.071;bmin=.144;bmax=.154;vmin=23.17;
        case 'Jabulani',a=.4839;vc=18.69;vs=1.377;bmin=.1413;bmax=.178;vmin=23.29;
        case 'Brazuca',a=.474;vc=12.92;vs=1;bmin=.1657;bmax=.2112;vmin=14.61;
        otherwise,error('未知实验球型');
    end
    d=(a-bmin)/(1+exp((U-vc)/vs))+bmin+ ...
        ((U-vmin)/(1+exp(-(U-vmin)/vs)))*(bmax-bmin)/(35-vmin);
end
qAir=.5*c.rho*pi*c.R^2*U^2;
if strcmp(c.liftModel,'teamgeist2018')
    d=d+.514*S;
    cl=1.15*S^.83*(.4927-min(d,.4927))/(.4927-(.144+.514*.19));
    if om>0,lift=qAir*cl*cross(omega/om,w/U);end
else
    % CL=k*S_perp with unit lift direction, written without a second sin(angle).
    lift=.5*c.rho*pi*c.R^3*c.liftSlope*cross(omega,w);
end
drag=-qAir*d*w/U;
end

function margin=domainMargin(t,s,c)
[~,~,U,S]=forces(t,s,c);
margin=1;
if strcmp(c.dragModel,'paper2018'),margin=min(U-7,35-U);end
if strcmp(c.liftModel,'teamgeist2018')
    if ~isfinite(S),margin=-1;else,margin=min([U-23.17,35-U,.4-S]);end
end
end

function [v,terminal,direction]=events(t,s,c)
ground=s(3)-c.R;if t==0&&ground==0,ground=eps;end
v=[s(1)-c.goalX;ground;domainMargin(t,s,c)];
terminal=[1;1;1];direction=[1;-1;-1];
end
