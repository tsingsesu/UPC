function family=penalty_family(obs,nominalCd,mode)
%PENALTY_FAMILY 网格剖面反演；保留全部候选/失败，不宣称全局包络或概率。
if nargin<3,mode='full';end
p=penalty_inverse('defaults');
if strcmp(mode,'quick'),winds=[-8 0 8];dens=[1.1 1.225 1.35];cds=unique([.25 nominalCd .5]);
else,winds=linspace(-8,8,7);dens=[1.1 1.225 1.35];cds=.15:.05:.6;end
n=numel(winds)*numel(dens)*numel(cds);
template=struct('windY',NaN,'density',NaN,'Cd',NaN,'theta',nan(1,6), ...
 'maxResidual',inf,'rmse',inf,'residual',nan(2,3),'exitflag',0,'boundaryHit',false, ...
 'accepted',false,'status','solver_failure','endpoint',nan(3,1),'tEnd',NaN, ...
 'initialSpeed',NaN,'kineticEnergy',NaN,'isOpening',false,'trajectory',zeros(0,3), ...
 'allStarts',[],'physicalAdmissible',false,'maxTransverseCoefficient',NaN,'rejectionReasons',{{}},'message','');
a=repmat(template,n,1);k=0;
for w=winds
 for rho=dens
  for cd=cds
   k=k+1;p.windY=w;p.density=rho;p.Cd=cd;
   a(k).windY=w;a(k).density=rho;a(k).Cd=cd;
   try
    f=penalty_inverse('fit',p,obs);s=penalty_inverse('propagate',f);
    a(k).theta=f.theta;a(k).maxResidual=f.fit.maxResidual;a(k).rmse=f.fit.rmse;
    a(k).residual=f.fit.residual;a(k).exitflag=f.fit.exitflag;a(k).boundaryHit=f.fit.boundaryHit;
    a(k).accepted=s.fit.feasible;a(k).physicalAdmissible=s.fit.physicalAdmissible;
    a(k).maxTransverseCoefficient=s.admissibility.maxTransverseCoefficient;a(k).rejectionReasons=s.fit.rejectionReasons;a(k).status=s.status;a(k).endpoint=s.endpoint(1:3);
    a(k).tEnd=s.tEnd;a(k).initialSpeed=s.initialSpeed;a(k).kineticEnergy=f.fit.kineticEnergy;
    a(k).isOpening=s.isOpening;a(k).allStarts=f.fit.allStarts;
    ts=unique([linspace(0,s.tEnd,100) .05 .1]);ts=ts(ts<=s.tEnd);
    vals=deval(s.sol,ts);a(k).trajectory=vals(1:3,:)';
   catch e
    a(k).message=e.message;
   end
  end
 end
 fprintf('[family] %d/%d profiles completed.\n',k,n);
end
summary=struct([]);
for tol=[.002 .005 .01 .02]
 ok=[a.exitflag]>0 & [a.maxResidual]<=tol & [a.physicalAdmissible];plane=ok&strcmp({a.status},'plane');
 pts=reshape([a(plane).endpoint],3,[]);speeds=[a(ok).initialSpeed];
 q=struct('tolerance',tol,'nAccepted',sum(ok),'nPlane',sum(plane), ...
  'nGround',sum(ok&strcmp({a.status},'ground')),'nTimeout',sum(ok&strcmp({a.status},'timeout')), ...
  'nOpening',sum(plane&[a.isOpening]),'yRange',span(pts(2,:)), ...
  'zRange',span(pts(3,:)),'speedRange',span(speeds));
 summary=[summary;q]; %#ok<AGROW>
end
accepted=[a.accepted]&strcmp({a.status},'plane');ids=find(accepted);pts=reshape([a(accepted).endpoint],3,[]);
best=NaN;pair=[];
for i=1:numel(ids)
 distances=sum((pts(2:3,:)-pts(2:3,i)).^2,1);[v,j]=max(distances);
 if isnan(best)||v>best,best=v;pair=[ids(i) ids(j)];end
end
family=struct('method','bounded least-squares profile representatives, componentwise residual screening', ...
 'scope','sampled inner subset, not exhaustive bounds or probability/confidence region', ...
 'defaultTolerance',.005,'windGrid',winds,'densityGrid',dens,'CdGrid',cds, ...
 'nCandidates',n,'nSolverFailures',sum([a.exitflag]<=0),'nBoundary',sum([a.boundaryHit]), ...
 'summary',summary,'physicalLimits',penalty_inverse('limits'),'nPhysicalRejected',sum(~[a.physicalAdmissible]&[a.exitflag]>0),'exampleIndices',pair,'exampleSeparation',sqrt(best),'members',a, ...
 'numericSettings',struct('fitRK4Step',.0025,'odeMaxStep',.002,'odeRelTol',1e-9, ...
 'odeAbsTol',1e-11,'fitStarts',2,'parameterBounds',[10 -20 0 -1.5 -1.5;40 20 25 1.5 1.5]));
sensitivity=struct([]);
for vmax=[35 40 45]
 for clmax=[.3 .4 .5 .6]
  ok=[a.exitflag]>0 & [a.maxResidual]<=.005 & [a.initialSpeed]<=vmax & [a.initialSpeed]>=10 & [a.maxTransverseCoefficient]<=clmax;
  sensitivity=[sensitivity;struct('speedCap',vmax,'coefficientCap',clmax,'accepted',sum(ok))]; %#ok<AGROW>
 end
end
family.admissibilitySensitivity=sensitivity;
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));img=fullfile(root,'img');
fig=figure('Visible','off','Color','w');ax=axes(fig);hold(ax,'on');
for i=ids
 r=a(i).trajectory;plot3(ax,r(:,1),r(:,2),r(:,3),'Color',[.75 .8 .86],'HandleVisibility','off');
end
colors=[.05 .3 .8;.85 .3 .05];
for j=1:numel(pair)
 r=a(pair(j)).trajectory;plot3(ax,r(:,1),r(:,2),r(:,3),'Color',colors(j,:),'LineWidth',1.6,'DisplayName',sprintf('Candidate %d',pair(j)));
end
plot3(ax,obs(:,1),obs(:,2),obs(:,3),'ko','MarkerFaceColor','k','DisplayName','Observations');
plot3(ax,11*[1 1 1 1],[-3.66 -3.66 3.66 3.66],[0 2.44 2.44 0],'k-','HandleVisibility','off');
xlabel(ax,'x (m)');ylabel(ax,'y (m)');zlabel(ax,'z (m)');grid(ax,'on');view(ax,-28,18);legend(ax,'Location','best');
exportgraphics(fig,fullfile(img,'全程约束轨迹簇.png'),'Resolution',200);close(fig);
fig=figure('Visible','off','Color','w');ax=axes(fig);hold(ax,'on');
scatter(ax,pts(2,:),pts(3,:),24,[a(accepted).Cd],'filled');cb=colorbar(ax);cb.Label.String='Drag coefficient C_D';
plot(ax,[-3.66 3.66 3.66 -3.66 -3.66],[0 0 2.44 2.44 0],'k-');
xlabel(ax,'y at x = 11 m (m)');ylabel(ax,'z at x = 11 m (m)');grid(ax,'on');axis(ax,'equal');
if ~isempty(pts)
 xlim(ax,[min(pts(2,:))-.15 max(pts(2,:))+.15]);ylim(ax,[min(pts(3,:))-.15 max(pts(3,:))+.15]);
else
 text(ax,0,1,'No compatible profile representative found','HorizontalAlignment','center');
end
exportgraphics(fig,fullfile(img,'全程约束球门散点.png'),'Resolution',200);close(fig);
T=table((1:n)',[a.windY]',[a.density]',[a.Cd]',[a.initialSpeed]',[a.kineticEnergy]', ...
 [a.maxResidual]',[a.accepted]',[a.exitflag]',{a.status}', ...
 'VariableNames',{'index','windY_mps','density_kgpm3','Cd','initialSpeed_mps','kineticEnergy_J', ...
 'maxResidual_m','accepted_5mm','exitflag','status'});
writetable(T,fullfile(root,'data','派生数据','全程约束参数簇.csv'));
end

function v=span(a)
if isempty(a),v=[NaN NaN];else,v=[min(a) max(a)];end
end
