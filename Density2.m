%function Density

clear all;
% http://github.com/rweigel/m-rsw/time
addpath('../../m-rsw/time')

% Read WGhourFS dataset from Kondrashov(2014)

% WGhourFS_72_13.txt from http:// on DATE
fnamem = 'data/WGhourFS_72_13.mat';
fnamet = 'data/WGhourFS_72_13.txt';

if exist(fnamem,'file')
    fprintf('Loading %s\n',fnamem)
    load('./data/WGhourFS_72_13.mat')
else

    if(~exist(fnamet,'file'))
        url = ['http://mag.gmu.edu/git-data/victoirv/Density/',fnamet];
        fprintf('Getting %s\n',url);
        urlwrite(url,fnamet);
    end

    fprintf('Reading %s\n',fnamet)

    KHead = textread(fnamet,'%s',28,'delimiter',',');
    KData = dlmread(fnamet,',',1,0);

    KVBs = 1/2*KData(:,6).*(abs(KData(:,5))-KData(:,5));
    KVBz = KData(:,6).*KData(:,5);
    KBs  = 1/2*(abs(KData(:,5))-KData(:,5));

    KDensity = KData(:,7);

    [Y,M,D] = ymd(KData(:,2),KData(:,1),'day');

    KTo = [Y(1),M(1),D(1),KData(1,3)];

    KDNo = datenum([Y(1),M(1),D(1)]); % Datenum of first row.
    KDN  = datenum([Y,M,D]);          % Datenum of each row.

    KH  = KData(:,3);

    fprintf('Saving %s\n',fnamem)
    save(fnamem,'KH','KDN','KVBs','KVBz','KBs','KDensity','KData','KHead')
end

%Get Denton data (mostly just for mass density)
fnamem = 'data/massdensity.mat';
fnamet = 'data/massdensity.txt';
if(exist(fnamem,'file'))
    fprintf('Loading %s\n',fnamem)
    load(fnamem)
else
    if(~exist(fnamet,'file'))
        url = ['http://mag.gmu.edu/git-data/victoirv/Density/',fnamet];
        fprintf('Getting %s\n',url);
        urlwrite(url,fnamet);
    end
    fprintf('Reading %s\n',fnamet);

    IN = dlmread(fnamet);
    IN(IN(:,1) ~= 7,:) = []; % Only use satellite 7 (6 should also work)
    IN(IN==9999) = NaN;    % Replace fill value with NaN.

    [Y,M,D] = ymd(IN(:,3),IN(:,2),'day');

    DTo = [Y(1),M(1),D(1),IN(1,4),IN(1,5)];

    DDNo = datenum([Y(1),M(1),D(1)]); % Datenum of first row.
    DDN  = datenum([Y,M,D]);          % Datenum of each row.

    DHM  = IN(:,4:5); % HR and Minute of day.

    DMN  = (DDN-DDNo)*24*60 + IN(:,4)*60 + IN(:,5); % Minute number relative to start

    [DMN,I] = unique(DMN);
    fprintf('Removed %d non-unique time values\n',length(DDN)-length(I))

    IN      = IN(I,:);
    DHM     = DHM(I,:);
    DDN     = DDN(I);

    % Column 88 is density, 10 is F10.7, 31 is BZ_sw, 32 is Vsw
    DF107 = IN(:,10);
    DBz   = IN(:,31);
    DVsw  = IN(:,32);
    DBs   = 1/2*IN(:,32).*(abs(DBz)-DBz);
    DMLT  = IN(:,8);
    DDensity = IN(:,88);

    DData = IN;
    fprintf('Saving %s\n',fnamem);
    save(fnamem,'DHM','DDN','DDensity','DF107','DBz','DMLT','DBs','DVsw','DData');
end

% Load F10.7 specifically 
fnamem = './data/omni2_1980_2000.mat';
fnamet = './data/omni2_1980_2000.lst';
if(exist(fnamem,'file'))
    fprintf('Loading %s\n',fnamem);   
    load(fnamem)
else
    if(~exist(fnamet,'file'))
        url = ['http://mag.gmu.edu/git-data/victoirv/Density/',fnamet];
        fprintf('Getting %s\n',url);
        urlwrite(url,fnamet);
    end

    OData = dlmread(fnamet);

    ODst  = OData(:,4);    
    OF107 = OData(:,5); 
    OF107(OF107 == 999.9) = NaN;

    [Y,M,D] = ymd(OData(:,2),OData(:,1),'day');
 
    OTo = [Y(1),M(1),D(1),OData(1,3)];

    ODNo = datenum([Y(1),M(1),D(1)]); % Datenum of first row.
    ODN  = datenum([Y,M,D]);          % Datenum of each row.

    OH  = OData(:,3);  % Hour relative to start.
    fprintf('Saving %s\n',fnamem);   
    save(fnamem,'OH','ODN','OF107','ODst','OData');
end

% Select overlapping windows
Io = max([ODN(1),DDN(1),KDN(1)]);
If = min([ODN(end),DDN(end),KDN(end)]);

IK = find(KDN>=Io & KDN <=If);
KDN = KDN(IK);
KH = KH(IK);
KVBs = KVBs(IK);
KVBz = KVBz(IK);
KBs  = KDensity(IK);
KData = KData(IK,:);
[KDN(1),KH(1);KDN(end),KH(end)]

ID = find(DDN>=Io & DDN <=If);
DDN = DDN(ID);
DHM = DHM(ID,:);
DDensity = DDensity(ID);
DF107 = DF107(ID);
DBz = DBz(ID);
DMLT = DMLT(ID);
DBs = DBs(ID);
DVsw = DVsw(ID);
[DDN(1),DHM(1,:);DDN(end),DHM(end,:)]

IO = find(ODN>=Io & ODN <=If);
ODN = ODN(IO);
OH = OH(IO);
OF107 = OF107(IO);
ODst = ODst(IO);
OData = OData(IO,:);
[ODN(1),OH(1);ODN(end),OH(end)]

n  = (DDN(end)-DDN(1)+1)*1440; % Number of minutes for fill array.

DDensityf       = NaN*ones(n,1); % Filled density on 1-min. grid.
Tf = [0:n-1];

To = (DDN-DDN(1))*1440 + DHM(:,1)*60 + DHM(:,2); % Minute since start day

DDensityf(To+1) = DDensity; % Replace NaNs w/ valid values.

% Reshape 1-min grid Density to have columns w/ 1-hr of data.
DDensityfr = reshape(DDensityf,60,length(DDensityf)/60);
DDensityMedian = nanmedian(DDensityfr,1);
DDensityNGood  = sum(~isnan(DDensityfr),1);
Tmedian = 60*[0:length(DDensityMedian)-1]; % Timestamps for median.

figure(1);clf;hold on;
    plot(To,DDensity,'b.','MarkerSize',30)
    plot(Tf,DDensityf,'g.','MarkerSize',20)
    plot(Tmedian,DDensityMedian,'k.','MarkerSize',10)
    for i = 1:100
        text(Tmedian(i),DDensityMedian(i)*1.01,num2str(DDensityNGood(i)));
    end
    legend('Original Data','Original on 1-min grid','Median over next 60 minutes')
    xlabel(['Minutes since ',datestr(DDN(1))])
    set(gca,'XLim',[0,1440-1])

I = find(ODst(1:end-1) > -30 & ODst(2:end) < -30);

k = 1;
for i = 1:length(I)
    a = I(i)-24;
    b = I(i)+24;
    if (a < 1),continue,end
    if (b > length(DDensityMedian)),break,end
    DDensityMedianStorm(k,:) = DDensityMedian(a:b);
    ODstStorm(k,:) = ODst(a:b);
    k = k+1;
end

t = [-24:24];
figure(2);clf;hold on;grid on;
    plot(t,nanmean(DDensityMedianStorm),'b')
    plot(t,nanmean(-ODstStorm),'g')
    xlabel('Time since onset [hrs]')
    legend('Density','-Dst')
