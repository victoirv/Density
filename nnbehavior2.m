function [testmean, testsd, net] = nnbehavior2(x,target,xtest,ytest,delays,loops)
if(nargin<4)
    delays=1;
    loops=10;
elseif(nargin<5)
    loops=10;
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
    
    for yloop=1:length(ytest)
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

