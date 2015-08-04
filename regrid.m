function [DNMedian,HMMedian,XMedian,XMean,NGood] = regrid(DN,HM,X)

% Number of minutes for fill array
n  = (DN(end)-DN(1)+1)*1440; 

% Minute since start day of available data values
To = (DN-DN(1))*1440 + HM(:,1)*60 + HM(:,2); 

% Arrays on 1-min. grid.
Xf  = NaN*ones(n,size(X,2));
Tf  = NaN*ones(n,1);
DNf = NaN*ones(n,1);
HMf = NaN*ones(n,2);

% Replace NaNs with valid values.
Tf(To+1)    = To;
DNf(To+1)   = DN;
HMf(To+1,:) = HM;
Xf(To+1,:)  = X;

% Place NaNs in all rows where last column of Xf (Density) is a NaN
If = find(isnan(Xf(:,end)));
Tf(If)    = NaN;
DNf(If)   = NaN;
HMf(If,:) = NaN;
for i = 1:size(Xf,2)-1
	Xf(If,i) = NaN;
end

% Reshape 1-min grid arrays have columns with 60 rows
DNfr = reshape(DNf,60,length(DNf)/60);
Tfr = reshape(Tf,60,length(Tf)/60);

DNMedian = nanmedian(DNfr,1);
TMedian  = nanmedian(Tfr,1);

for i = 1:2
	HMfr{i} = reshape(HMf(:,i),60,length(HMf(:,i))/60);
	% Find median and mean of non-nan points in 60-minute windows.
	HMMedian(:,i) = nanmedian(HMfr{i},1);
end

Nc = size(Xf,2);
for i = 1:Nc
	Xfr{i}  = reshape(Xf(:,i),60,length(Xf(:,i))/60);

	if (i == Nc)
		% Number of valid points used in means and medians
		% is equal to number of non-NaNs in last column (actually all
		% because of code above).
		NGood = sum(~isnan(Xfr{Nc}),1);
	end

	% Find median and mean of non-nan points in 60-minute windows.
	XMedian(:,i) = nanmedian(Xfr{i},1);
	XMean(:,i)   = nanmean(Xfr{i},1);
end
if (0)
figure(1);clf;hold on;
    plot(To,X(:,end),'b.','MarkerSize',30)
    plot(Tf,Xf(:,end),'g.','MarkerSize',20)
    plot(TMedian,XMedian(:,end),'k.','MarkerSize',10)
    for i = 1:100
        text(TMedian(i),XMedian(i,end)*1.01,num2str(NGood(i)));
    end
    legend('Original Data','Original on 1-min grid','Median in 60 minute window')
    for i = 600:60:1400
    	plot([i,i],[3,11],'k')
    end
    xlabel(['Minutes since ',datestr(DN(1))])
    set(gca,'XLim',[0,1440-1])
end