function output=NNBinaryOnset(x,target,figname)
if(nargin<3)
    figname='';
end


delays=1;

if(length(size(x))==3) %If sending in multiple variables for each storm and time, just add them as a vector of variables for each storm
    x=reshape(x,length(x),[]);
end

Z=[x target];
for col=1:min(size(Z))
    Z(isnan(Z(:,col)),:)=[];
end
x=Z(:,1:end-1);
y=Z(:,end);




inputDelays = 1:delays;
    feedbackDelays = 1:delays;
    hiddenLayerSize = 10;
    
    net=patternnet(hiddenLayerSize);
    net=train(net,x',y'); %Look at confusion matrix. Not great.
    
    output=net(x');
    
    
    plotconfusion(y',output)
    if(~isempty(figname))
    print('-depsc2',sprintf('figures/NNBinaryOnset-%s.eps',figname)); 
    print('-dpng','-r200',sprintf('figures/PNGs/NNBinaryOnset-%s.png',figname)); 
    end
    
    %{
    net = narxnet(inputDelays,feedbackDelays,hiddenLayerSize);
    net = removedelay(net);
    [inputs,inputStates,layerStates,targets] = preparets(net,tonndata(x,false,false),{},tonndata(y,false,false));
    
    net.divideFcn = 'divideblock';
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    

    
    [net,tr,Ys Es Xf Af] = train(net,inputs,targets,inputStates,layerStates);
    
    outputs = net(inputs,inputStates,layerStates);
    %}