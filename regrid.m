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

% Replace NaNs with valid values. f = filled.
Tf(To+1)    = To;
DNf(To+1)   = DN;
HMf(To+1,:) = HM;
Xf(To+1,:)  = X;

% Reshape 1-min grid arrays have columns with 60 rows. r = reshaped.
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