function [testmean, testsd] = nnbehavior(x,target,xtest,delays,loops,plotvars,AXmap,satnum)
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
if(diff(size(xtest))>0), xtest=xtest'; end

if(length(size(x))==3) %If sending in multiple variables for each storm and time, just add them as a vector of variables for each storm
    x=reshape(x,length(x),[]);
end

Z=[x target];
for col=1:min(size(Z))
    Z(isnan(Z(:,col)),:)=[];
end
x=Z(:,1:end-1);
target=Z(:,end);

xtest=linspace(nanmean(x)-nanstd(x),nanmean(x)+nanstd(x))'; %Only test across points still valid after culling

hiddenLayerSize = 10;


net=timedelaynet(1:delays,hiddenLayerSize);
[testinputs,testinputStates,testlayerStates] = preparets(net,tonndata(xtest,false,false));
outputs=zeros(loops,length(xtest)-1);

parfor loop=1:loops
    net=timedelaynet(1:delays,hiddenLayerSize);
    [inputs,inputStates,layerStates,targets] = preparets(net,tonndata(x,false,false),tonndata(target,false,false));
    
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;  
    net.divideParam.testRatio = 15/100;
    net.trainParam.showWindow=0;

    net = train(net,inputs,targets,inputStates,layerStates);   
    
    outputs(loop,:) = fromnndata(net(testinputs,testinputStates,testlayerStates),1,1,0);
end

if(loops>1)
    testmean=squeeze(nanmean(outputs,1));
    testsd=squeeze(nanstd(outputs,1));
else
    testmean=squeeze(outputs);
    testsd=nan(size(testmean));
end

if(plotv)
    figure; plot(xtest(2:end),testmean,'.')
    hold on; plot(xtest(2:end),[testmean+testsd; testmean-testsd],'r.')
    xlabel(plotvars{1})
    ylabel(plotvars{2})
   % axis([AXmap(plotvars{1}) AXmap(plotvars{2})])
    title(sprintf('Mean nonlinear predicted %s over %d loops for GOES %d',plotvars{2},loops,satnum))
        filename=sprintf('figures/NN%s-GOES%d.',strjoin(plotvars,'-'),satnum);
    filename=regexprep(filename,'[^a-zA-Z0-9/\-]','');
    print('-depsc2',strcat(filename,'.eps'))
    print('-dpng',strcat(filename,'.png'))
end

