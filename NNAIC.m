function NNAIC(x,target,delays)


Z=[x target];
for col=1:min(size(Z))
    Z(isnan(Z(:,col)),:)=[];
end
x=Z(:,1:end-1);
target=Z(:,end);


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