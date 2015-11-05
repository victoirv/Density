function CCM = nntest(x,y,delays,silent)
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

Z=[x y];
for col=1:min(size(Z))
    Z(isnan(Z(:,col)),:)=[];
end
x=Z(:,1:end-1);
y=Z(:,end);

inputDelays = 1:delays;
feedbackDelays = 1:delays;
hiddenLayerSize = 10;
net = narxnet(inputDelays,feedbackDelays,hiddenLayerSize);
net = removedelay(net);
[inputs,inputStates,layerStates,targets] = preparets(net,tonndata(x,false,false),{},tonndata(y,false,false));

net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

if(silent)
    net.trainParam.showWindow=0;
end

[net,tr,Ys Es Xf Af] = train(net,inputs,targets,inputStates,layerStates);

outputs = net(inputs,inputStates,layerStates);
errors = gsubtract(targets,outputs);
performance = perform(net,targets,outputs);

%since we're about to shift the output by 1 to predict one step ahead, cut
%off last index if it's the last data point. "targets" variable is already
%shifted by 1 via removedelay() so it coincides for comparison
if(max(tr.trainInd)>=length(y)), tr.trainInd(end)=[]; end 
if(max(tr.testInd)>=length(y)), tr.testInd(end)=[]; end
if(max(tr.valInd)>=length(y)), tr.valInd(end)=[]; end

Ztr=[x(tr.trainInd,:) y(tr.trainInd+1)];
%{
for col=1:min(size(Ztr))
    Ztr(isnan(Ztr(:,col)),:)=[];
end
%}
coef=Ztr(:,1:end-1)\Ztr(:,end);
cctr=corrcoef(Ztr(:,1:end-1)*coef,Ztr(:,end)); cctr=cctr(1,2);
nntr=corrcoef(fromnndata(targets(tr.trainInd),1,0,0),fromnndata(outputs(tr.trainInd),1,0,0),'rows','pairwise'); nntr=nntr(1,2);

Zt=[x(tr.testInd,:) y(tr.testInd+1)];
for col=1:min(size(Zt))
    Zt(isnan(Zt(:,col)),:)=[];
end
coef=Zt(:,1:end-1)\Zt(:,end);
cct=corrcoef(Zt(:,1:end-1)*coef,Zt(:,end)); cct=cct(1,2);
nnt=corrcoef(fromnndata(targets(tr.testInd),1,0,0),fromnndata(outputs(tr.testInd),1,0,0),'rows','pairwise'); nnt=nnt(1,2);

Zv=[x(tr.valInd,:) y(tr.valInd+1)];
for col=1:min(size(Zv))
    Zv(isnan(Zv(:,col)),:)=[];
end
coef=Zv(:,1:end-1)\Zv(:,end);
ccv=corrcoef(Zv(:,1:end-1)*coef,Zv(:,end)); ccv=ccv(1,2);
nnv=corrcoef(fromnndata(targets(tr.valInd),1,0,0),fromnndata(outputs(tr.valInd),1,0,0),'rows','pairwise'); nnv=nnv(1,2);

CCM=[cctr cct ccv; nntr nnt nnv]; %Correlation matrix for plotting later

if(~silent)
fprintf(' CC(#valid/tot)      \t Linear\t Neural \n Train(%d/%d) \t\t %2.3f \t %2.3f \n Test(%d/%d) \t\t %2.3f \t %2.3f \n Validate(%d/%d) \t %2.3f \t %2.3f \n',...
    sum(~isnan(fromnndata(targets(tr.trainInd),1,0,0))),length(tr.trainInd),cctr, nntr, sum(~isnan(fromnndata(targets(tr.testInd),1,0,0))),length(tr.testInd), cct,nnt,sum(~isnan(fromnndata(targets(tr.valInd),1,0,0))),length(tr.valInd), ccv,nnv)
end

