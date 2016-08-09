function [testmean, testsd] = nnbehavior2(x,target,delays,loops,plotvars,AXmap,satnum)
if(nargin<4)
    delays=1;
    loops=10;
elseif(nargin<5)
    loops=10;
end

plotv=0;
if(exist('plotvars','var'))
    plotv=1;
end

UnitMap=containers.Map({'B_z','V_{sw}','F_{10.7}','\rho_{eq}','DoY','\rho_{sw}','D_{st}','MLT'},{'nT','m/s','s.f.u.','amu/cm^3','days','amu/cm^3','nT','hour'});


%Make column vectors
if(diff(size(x))>0), x=x'; end

if(length(size(x))==3) %If sending in multiple variables for each storm and time, just add them as a vector of variables for each storm
    x=reshape(x,length(x),[]);
end

Z=[x target];
for col=1:min(size(Z))
    Z(isnan(Z(:,col)),:)=[];
end
x=Z(:,1:end-1);
target=Z(:,end);

xtest=linspace(min(x(:,1)),max(x(:,1)))'; %Only test across points still valid after culling
ytest=linspace(min(x(:,2)),max(x(:,2)))';

hiddenLayerSize = 10;
outputs=zeros(loops,length(ytest),length(xtest)-1);

for loop=1:loops
    net=timedelaynet(1:delays,hiddenLayerSize);
    [inputs,inputStates,layerStates,targets] = preparets(net,tonndata(x,false,false),tonndata(target,false,false));
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;  
    net.divideParam.testRatio = 15/100;
    net.trainParam.showWindow=0;
    [net, tr] = train(net,inputs,targets,inputStates,layerStates); 
    
    coefs(:,loop)=[x(tr.trainInd,:) ones(length(tr.trainInd),1)]\target(tr.trainInd);
    
    parfor yloop=1:length(ytest)
        [testinputs,testinputStates,testlayerStates] = preparets(net,tonndata([xtest ones(length(xtest),1).*ytest(yloop)],0,0) ); 
        outputs(loop,yloop,:) = fromnndata(net(testinputs,testinputStates,testlayerStates),1,1,0);
    end
    
end

if(loops>1)
    testmean=squeeze(nanmean(outputs,1));
    testsd=squeeze(nanstd(outputs,1));
else
    testmean=squeeze(outputs);
    testsd=nan(size(testmean));
end

if(plotv)
    scattersize=round((length(xtest)*length(ytest))*0.005); %So no point takes up more than ~0.1% of the graph
    figure; surf(xtest(2:end),ytest,testmean,'EdgeColor','none','LineStyle','none','FaceLighting','phong') %Even though phong is deprecated, it's the only one that plots without corruption
    
    ylabel(sprintf('%s (%s)',plotvars{2},UnitMap(plotvars{2})))
    xlabel(sprintf('%s (%s)',plotvars{1},UnitMap(plotvars{1})))
    %axis([AXmap(plotvars{1}) AXmap(plotvars{2}) AXmap(plotvars{3})])
    %caxis(AXmap('\rho_{eq}'))    
    view(0,90)
    title(sprintf('Mean nonlinear predicted %s over %d loops with GOES %d',plotvars{3},loops,satnum))
    colorbar
    hold on; scatter3(x(:,1),x(:,2),repmat(500,1,length(x(:,1))),(target-min(target(:)))./max(target(:)).*scattersize+1,'k')
    %set(gca,'YTick',linspace(yr(1),yr(2),10))
    
    ax=[xlim; ylim; zlim; caxis];
    filename=sprintf('figures/CorrPlots/NN%s-GOES%d.',strjoin(plotvars,'-'),satnum);
    filename=regexprep(filename,'[^a-zA-Z0-9/\-]','');
    print('-depsc2','-zbuffer','-r300',strcat(filename,'.eps'))
    print('-dpng',strcat(filename,'.png'))
    
    
    
    figure; surf(xtest(2:end),ytest,testsd,'EdgeColor','none','LineStyle','none','FaceLighting','phong') %Even though phong is deprecated, it's the only one that plots without corruption
    ylabel(sprintf('%s (%s)',plotvars{2},UnitMap(plotvars{2})))
    xlabel(sprintf('%s (%s)',plotvars{1},UnitMap(plotvars{1})))
    %axis([AXmap(plotvars{1}) AXmap(plotvars{2})])
    title(sprintf('Standard deviation of %d loops predicting %s with GOES %d',loops,plotvars{3},satnum))
    colorbar
    axis(reshape(ax(1:3,:)',1,[]))
    view(0,90)
        filename=sprintf('figures/CorrPlots/NN%s-sd-GOES%d.',strjoin(plotvars,'-'),satnum);
    filename=regexprep(filename,'[^a-zA-Z0-9/\-]','');
    print('-depsc2','-zbuffer','-r300',strcat(filename,'.eps'))
    print('-dpng',strcat(filename,'.png'))
    
    
    coef=[x ones(length(x),1)]\target;
    coef=nanmedian(coefs,2);
    [Zx, Zy]=meshgrid(xtest,ytest);
    Z=Zx.*coef(1)+Zy.*coef(2)+1.*coef(3);
    cc=corrcoef([x ones(length(x),1)]*coef,target);
    
    
    
    figure;
    surf(xtest,ytest,Z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
    
    ylabel(sprintf('%s (%s)',plotvars{2},UnitMap(plotvars{2})))
    xlabel(sprintf('%s (%s)',plotvars{1},UnitMap(plotvars{1})))
    %axis([AXmap(plotvars{1}) AXmap(plotvars{2}) AXmap(plotvars{3})])
   
    %caxis(AXmap('\rho_{eq}'))
    view(0,90)
    grid on
    %set(gca,'YTick',linspace(yr(1),yr(2),10))
    %set(gca,'XTick',linspace(xr(1),xr(2),10))
    title(sprintf('Mean linear predicted %s over %d loops with GOES %d - CC%2.2f',plotvars{3},loops,satnum,cc(1,2)))
    colorbar
    hold on; scatter3(x(:,1),x(:,2),repmat(500,1,length(x(:,1))),(target-min(target(:)))./max(target(:)).*scattersize+1,'k')
     axis(reshape(ax',1,[]))
    filename=sprintf('figures/CorrPlots/Linear%s-GOES%d.',strjoin(plotvars,'-'),satnum);
    filename=regexprep(filename,'[^a-zA-Z0-9/\-]','');
    print('-depsc2','-zbuffer','-r300',strcat(filename,'.eps'))
    print('-dpng',strcat(filename,'.png'))
    
end

