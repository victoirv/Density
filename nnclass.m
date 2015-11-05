function nnclass(x,y,delays,silent)
if(nargin<2)
    x=rand(1,1000);
    y=[0.4.*x(2:end)-0.2.*x(1:end-1) 0];
    delays=1;
    silent=0;
elseif(nargin<3)
    delays=1;
    silent=0;
elseif(nargin<4)
    silent=0;
end
if(length(size(x))==3) %If sending in multiple variables for each storm and time, just add them as a vector of variables for each storm
    x=reshape(x,length(x),[]);
end

inputDelays = 1:delays;
feedbackDelays = 1:delays;
hiddenLayerSize = 10;

%Make it a classification problem
split=quantile(y,[0.33 0.66]);
yc=y.*0;
yc(y>split(2))=1; yc(y<split(1))=-1;
netc = narxnet(inputDelays,feedbackDelays,hiddenLayerSize);
netc.divideParam.trainRatio = 70/100;
netc.divideParam.valRatio = 15/100;
netc.divideParam.testRatio = 15/100;
netc = removedelay(netc);
[xc,xic,aic,tc] = preparets(netc,tonndata(x,false,false),{},tonndata(yc,false,false));

[netc,trc] = train(netc,xc,tc,xic,aic);
[Yc,Xf,Af] = netc(xc,xic,aic);
errorss=gsubtract(tc,Yc);
performancec = perform(netc,tc,Yc);
