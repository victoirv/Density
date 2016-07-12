function output=NNBinaryOnset(x,target,figname,labels)
if(nargin<3)
    figname='';
end
if(nargin<4)
    labels='';
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
net=train(net,x',y',[],[],y'.*log(length(y)/sum(y))+1); %Look at confusion matrix. Not great.

output=net(x');


plotconfusion(y',output)
if(~isempty(figname))
    print('-depsc2',sprintf('figures/NNBinaryOnset-%s.eps',figname));
    print('-dpng','-r200',sprintf('figures/PNGs/NNBinaryOnset-%s.png',figname));
end



%Make histograms if labels are passed in
if(~isempty(labels))
    
    if(findstr('full',figname))
        for i=1:length(labels)
        x2(:,i)=nanmedian(x(:,i:length(labels):end),2);
        end
        x=x2;
    end
        [bina,binax]=hist(x(y & 1,1)); %Not sure why I need the & 1, but it doesn't seem to think target is 1s and 0s otherwise
        [bin,binx]=hist(x(round(output) & y',1));
        [bin2,binx2]=hist(x(~round(output) & y',1));
        [bin3,binx3]=hist(x(round(output) & ~y',1));


        if(min(size(x))==1)
            figure; stairs(binax,bina,'k','LineWidth',1.2);
            hold on; stairs(binx,bin); stairs(binx2,bin2,'r'); stairs(binx3,bin3,'m');
            legend('Actual','True Positive','False Negative','False Positive')
            xlabel(labels)
        else
            figure; subplot(min(size(x)),1,1)
            stairs(binax,bina,'k','LineWidth',1.2);
            hold on; stairs(binx,bin); stairs(binx2,bin2,'r'); stairs(binx3,bin3,'m');
            legend('Actual','True Positive','False Negative','False Positive')
            xlabel(labels{1})

            for i=2:min(size(x))
                [bina,binax]=hist(x(y & 1,i));
                [bin,binx]=hist(x(round(output) & y',i));
                [bin2,binx2]=hist(x(~round(output) & y',i));
                [bin3,binx3]=hist(x(round(output) & ~y',i));
                subplot(min(size(x)),1,i)
                stairs(binax,bina,'k','LineWidth',1.2);
                hold on; stairs(binx,bin); stairs(binx2,bin2,'r'); stairs(binx3,bin3,'m');
                legend('Actual','True Positive','False Negative','False Positive')
                xlabel(labels{i})
            end

        end
    %end
    
    if(~isempty(figname))
        orient tall;
        print('-depsc2',sprintf('figures/NNBinaryOnset-%s-hist.eps',figname));
        print('-dpng','-r200',sprintf('figures/PNGs/NNBinaryOnset-%s-hist.png',figname));
    end
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