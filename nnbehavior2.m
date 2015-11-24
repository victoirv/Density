function [testmean, testsd] = nnbehavior2(x,target,delays,loops,plotvars,satnum)
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
    net = train(net,inputs,targets,inputStates,layerStates); 
    
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
    figure; surf(xtest(2:end),ytest,testmean,'EdgeColor','none','LineStyle','none','FaceLighting','phong') %Even though phong is deprecated, it's the only one that plots without corruption
    view(0,90)
    ylabel(plotvars{2})
    xlabel(plotvars{1})
    title(sprintf('Mean predicted %s over %d loops with GOES %d',plotvars{3},loops,satnum))
    colorbar
    hold on; scatter3(x(:,1),x(:,2),repmat(100000,1,length(x(:,1))),target,'k')
    print('-dpng',sprintf('figures/NN%s-GOES%d.png',strjoin(plotvars,'-'),satnum))
    
    figure; surf(xtest(2:end),ytest,testsd,'EdgeColor','none','LineStyle','none','FaceLighting','phong') %Even though phong is deprecated, it's the only one that plots without corruption
    view(0,90)
    ylabel(plotvars{2})
    xlabel(plotvars{1})
    title(sprintf('Standard deviation of %d loops predicting %s with GOES %d',loops,plotvars{3},satnum))
    colorbar
    print('-dpng',sprintf('figures/NN%s-sd-GOES%d.png',strjoin(plotvars,'-'),satnum))
    
    coef=[x ones(length(x),1)]\target;
    [Zx, Zy]=meshgrid(xtest,ytest);
    Z=Zx.*coef(1)+Zy.*coef(2)+1.*coef(3);
    figure;
    surf(xtest,ytest,Z)
    view(0,90)
    ylabel(plotvars{2})
    xlabel(plotvars{1})
    xlim([-10,10])
    zlim([0,50])
    title(sprintf('Linear predicted %s over %d loops with GOES %d',plotvars{3},loops,satnum))
    colorbar
    hold on; scatter3(x(:,1),x(:,2),repmat(50,1,length(x(:,1))),target,'k')
    print('-dpng',sprintf('figures/Linear%s-GOES%d.png',strjoin(plotvars,'-'),satnum))
    
end

