%function Density

clear all;
% http://github.com/rweigel/m-rsw/time
if ~exist('m-rsw')
    system('git clone http://github.com/rweigel/m-rsw')
    addpath('./m-rsw/time')
end
if ~exist('data')
    mkdir('data');
end

%for sat = [2,3,5,6,7]
for sat = [6]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read WGhourFS dataset from Kondrashov (2014)
    % WGhourFS_72_13.txt from
    % http://research.atmos.ucla.edu/tcd//dkondras/SolarWind.html
    fnamem = 'data/WGhourFS_72_13.mat';
    fnamet = 'data/WGhourFS_72_13.txt';

    if exist(fnamem,'file')
        %fprintf('Loading %s\n',fnamem)
        load(fnamem)
    else
        if(~exist(fnamet,'file'))
            url = ['http://mag.gmu.edu/git-data/victoirv/Density/',fnamet];
            fprintf('Getting %s\n',url);
            urlwrite(url,fnamet);
            fprintf('Wrote %s\n',fnamet);
        end

        fprintf('Reading %s\n',fnamet)

        KHead = textread(fnamet,'%s',28,'delimiter',',');
        KData = dlmread(fnamet,',',1,0);

        KBz  = KData(:,5);
        KVsw = KData(:,6);
        KDensitySW = KData(:,7);
        KDst = KData(:,15);

        [Y,M,D] = ymd(KData(:,2),KData(:,1),'day');

        KTo = [Y(1),M(1),D(1),KData(1,3)];

        KDNo = datenum([Y(1),M(1),D(1)]); % Datenum of first row.
        KDN  = datenum([Y,M,D]);          % Datenum of each row.

        KH  = KData(:,3);

        fprintf('Saving %s\n',fnamem)
        save(fnamem,'KH','KDN','KBz','KVsw','KDensitySW','KDst','KData','KHead')
    end


    % Load data from OMNI2 dataset
    % http://omniweb.gsfc.nasa.gov/form/dx1.html
    % Accessed on 07/29/2015

    %fnamem = './data/omni2_1980_2000.mat';
    %fnamet = './data/omni2_1980_2000.lst';
    fnamem = 'data/omni2_6477.mat';
    fnamet = 'data/omni2_6477.lst';
    if(exist(fnamem,'file'))
        %fprintf('Loading %s\n',fnamem);   
        load(fnamem);
    else
        if(~exist(fnamet,'file'))
            url = ['http://mag.gmu.edu/git-data/victoirv/Density/',fnamet];
            fprintf('Getting %s\n',url);
            urlwrite(url,fnamet);
            fprintf('Wrote %s\n',fnamet);
        end

        OData = dlmread(fnamet);

        OBz = OData(:,4); 
        OBz(OBz == 999.9) = NaN; 
        ODensitySW = OData(:,5); 
        ODensitySW(ODensitySW == 999.9) = NaN; 
        OVsw = OData(:,6); 
        OVsw(OVsw == 9999) = NaN; 
        OKp = OData(:,7);
        OKp(OKp == 99) = NaN; 
        ODst  = OData(:,8);
        ODst(ODst == 99999) = NaN;
        OF107 = OData(:,9); 
        OF107(OF107 == 999.9) = NaN;

        [Y,M,D] = ymd(OData(:,2),OData(:,1),'day');
     
        OTo = [Y(1),M(1),D(1),OData(1,3)];

        ODNo = datenum([Y(1),M(1),D(1)]); % Datenum of first row.
        ODN  = datenum([Y,M,D]);          % Datenum of each row.

        OH  = OData(:,3);  % Hour relative to start.
        fprintf('Saving %s\n',fnamem);
        save(fnamem,'OH','ODN','OBz','OVsw','ODensitySW','ODst','OF107','OData');
    end

    fprintf('Satellite: GOES-%d\n',sat);

    % Denton data
    % http://www.dartmouth.edu/~rdenton/Data/DentonTakahashiGOES1980-1991MassDensityWithHeader.txt
    % Has garbage characters starting on 25100.

    fnamem = sprintf('data/massdensity_sat_%d.mat',sat);
    fnamet = 'data/massdensity.txt';
    if(exist(fnamem,'file'))
        %fprintf('Loading %s\n',fnamem)
        load(fnamem)
    else
        if(~exist(fnamet,'file'))
            url = ['http://mag.gmu.edu/git-data/victoirv/Density/',fnamet];
            fprintf('Getting %s\n',url);
            urlwrite(url,fnamet);
            fprintf('Wrote %s\n',fnamet);
        end
        fprintf('Reading %s\n',fnamet);
     
        IN = dlmread(fnamet);
        fprintf('Unique satellite numbers: ');
        fprintf('%d ',unique(IN(:,1)))
        IN(IN(:,1) ~= sat,:) = []; % Only use satellite sat
        IN(IN==9999) = NaN;      % Replace fill value with NaN.

        [Y,M,D] = ymd(IN(:,3),IN(:,2),'day');

        DTo = [Y(1),M(1),D(1),IN(1,4),IN(1,5)];

        DDNo = datenum([Y(1),M(1),D(1)]); % Datenum of first row.
        DDN  = datenum([Y,M,D]);          % Datenum of each row.

        DHM  = IN(:,4:5); % HR and Minute of day.

        % Minute number relative to start
        DMN  = (DDN-DDNo)*24*60 + IN(:,4)*60 + IN(:,5); 

        [DMN,I] = unique(DMN);
        fprintf('Removed %d non-unique time values\n',length(DDN)-length(I))

        IN      = IN(I,:);
        DHM     = DHM(I,:);
        DDN     = DDN(I);

        % Column 88 is mass density, 10 is F10.7, 31 is Bz, 32 is Vsw
        DMLT  = IN(:,8);
        DF107 = IN(:,10);
        DDst  = IN(:,13);
        DBz   = IN(:,31);
        DVsw  = IN(:,32);
        DDensitySW = IN(:,33);
        DDensity = IN(:,88);

        DData = IN;
        fprintf('Saving %s\n',fnamem);
        save(fnamem,'DHM','DDN','DBz','DVsw','DDensitySW',...
                'DDst','DF107','DDensity','DMLT','DData');
    end

    % Select overlapping windows
    Io = max([ODN(1),DDN(1),KDN(1)]);
    Io = datenum('1989-01-01');
    If = min([ODN(end),DDN(end),KDN(end)]);

    fprintf('Overlapping time window is %s through %s\n',datestr(Io),datestr(If));

    IK   = find(KDN>=Io & KDN <=If);
    KDN  = KDN(IK);
    KH   = KH(IK);
    KBz  = KBz(IK);
    KVsw = KVsw(IK);
    KDensitySW = KDensitySW(IK);
    KDst  = KDst(IK);
    KData = KData(IK,:);

    ID  = find(DDN>=Io & DDN <=If);
    DDN = DDN(ID);
    DHM = DHM(ID,:);
    DBz = DBz(ID);
    DVsw = DVsw(ID);
    DDensitySW = DDensitySW(ID);
    DF107 = DF107(ID);
    DDst  = DDst(ID);
    DDensity = DDensity(ID);
    DMLT = DMLT(ID);
    %[DDN(1),DHM(1),DF107(1);DDN(end),DHM(end,:),DF107(end)]

    fprintf('Computing means and medians in 1-hour windows for Denton data.\n')
    [DDNMedian,DHMMedian,DXMedian,DXMean,DNGood] = ....
        regrid(DDN,DHM,[DBz,DVsw,DDensitySW,DDst,DF107,DMLT,DDensity]);

    DBzMedian = DXMedian(:,1);
    DVswMedian = DXMedian(:,2);
    DDensitySWMedian = DXMedian(:,3);
    DDstMedian = DXMedian(:,4);
    DF107Median = DXMedian(:,5);
    DMLTMedian = DXMedian(:,6);
    DDensityMedian = DXMedian(:,7);

    DBzMean = DXMean(:,1);
    DVswMean = DXMean(:,2);
    DDensitySWMean = DXMean(:,3);
    DDstMean = DXMean(:,4);
    DF107Mean = DXMean(:,5);
    DMLTMean = DXMean(:,6);
    DDensityMean = DXMean(:,7);

    IO = find(ODN>=Io & ODN <=If);
    ODN = ODN(IO);
    OH = OH(IO);
    OBz = OBz(IO);
    OVsw = OVsw(IO);
    OF107 = OF107(IO);
    ODst = ODst(IO);
    OData = OData(IO,:);
    %[ODN(1),OH(1);ODN(end),OH(end)]

    tmp = corrcoef(KBz,OBz,'rows','complete');
    fprintf('cc(KBz,OBz)                    = %.3f\n',tmp(2));
    tmp = corrcoef(DBzMedian,OBz,'rows','complete');
    fprintf('cc(DBzMedian,OBz)              = %.3f\n',tmp(2));
    tmp = corrcoef(DBzMean,OBz,'rows','complete');
    fprintf('cc(DBzMean,OBz)                = %.3f\n',tmp(2));

    tmp = corrcoef(KVsw,OVsw,'rows','complete');
    fprintf('cc(KVsw,OVsw)                  = %.3f\n',tmp(2));
    tmp = corrcoef(DVswMedian,OVsw,'rows','complete');
    fprintf('cc(DVswMedian,OVsw)            = %.3f\n',tmp(2));
    tmp = corrcoef(DVswMean,OVsw,'rows','complete');
    fprintf('cc(DVswMean,OVsw)              = %.3f\n',tmp(2));

    tmp = corrcoef(OF107,DF107Median,'rows','complete');
    fprintf('cc(OF107,DF107Median)          = %.3f\n',tmp(2));
    tmp = corrcoef(OF107,DF107Mean,'rows','complete');
    fprintf('cc(OF107,DF107Mean)            = %.3f\n',tmp(2));

    tmp = corrcoef(KDst,ODst,'rows','complete');
    fprintf('cc(KDst,ODst)                  = %.3f\n',tmp(2));
    tmp = corrcoef(ODst,DDstMedian,'rows','complete');
    fprintf('cc(ODst,DDstMedian)            = %.3f\n',tmp(2));
    tmp = corrcoef(ODst,DDstMean,'rows','complete');
    fprintf('cc(ODst,DDstMean)              = %.3f\n',tmp(2));

    tmp = corrcoef(OF107,DDensityMedian,'rows','complete');
    fprintf('cc(OF107,DDensityMedian)       = %.3f\n',tmp(2));
    tmp = corrcoef(OF107,DDensityMean,'rows','complete');
    fprintf('cc(OF107,DDensityMean)         = %.3f\n',tmp(2));

    tmp = corrcoef(OF107,log10(DDensityMedian),'rows','complete');
    fprintf('cc(OF107,log10(DDensityMedian) = %.3f\n',tmp(2));
    tmp = corrcoef(OF107,log10(DDensityMean),'rows','complete');
    fprintf('cc(OF107,log10(DDensityMean)   = %.3f\n',tmp(2));

    %I = find(DDstMedian(1:end-1) > -50 & DDstMedian(2:end) < -50);
    %I = find(ODst(1:end-1) > -50 & ODst(2:end) < -50);
    I = find(KDst(1:end-1) > -50 & KDst(2:end) <= -50);
    %I = find(KDst(1:end-1) > -50 & KDst(2:end) < -50);
    clear *Storm*

    Nd = 2;
    k = 1;
    for i = 1:length(I)-1
        a = I(i)-24*Nd+12;
        b = I(i)+24*Nd-12;
        %if ( DMLT(I(i)) < 6 || DMLT(I(i)) > 12 ),continue,end
        %if (I(i+1) < I(i) + 12),continue,end
        if (a < 1),continue,end
        if (b > length(DDensityMedian)),break,end
        fprintf('Event at %s %.1f %.1f %.1f\n',...
            datestr(Io+I(i)/24),KDst(I(i)),KDst(I(i)+1),KDst(I(i)+2))
        DDensityStorm(k,:) = DDensityMedian(a:b);
        DDstStorm(k,:)  = DDstMedian(a:b);
        ODstStorm(k,:)  = ODst(a:b);
        DF107Storm(k,:)  = DF107Median(a:b);
        DNGoodStorm(k,:) = DNGood(a:b);
        k = k+1;
    end

    t = [-24*Nd+12:24*Nd-12];

    for j = 1:size(DDensityStorm,2)/24 % Time in event
        a = 1+24*(j-1);
        b = 24*j;
        tmp1 = DDensityStorm(:,[a:b]);
        tmp2 = nanmedian(tmp1,1);
        tmp2n = sum(~isnan(tmp1(:,1)));
        DDensityStorm1Day1(j) = nanmedian(tmp1(:));
        DDensityStorm1Day2(j) = nanmedian(tmp2);
        N1D1(j) = sum(~isnan(tmp1(:)));        
        N1D2(j) = tmp2n;     

        %DDensityStorm1Day(j) = nanmean(tmp(:));
        tc(j) = mean(t(a:b));
    end

    if (length(I) > 0)
        figure(1);clf;hold on;grid on;
            plot(t,nanmedian(DDensityStorm),'b','LineWidth',3)
            plot(t,nanmedian(-DDstStorm),'g','LineWidth',3)
            plot(t,nanmean(-ODstStorm),'k','LineWidth',3)
            plot(t,nanmedian(DF107Storm)/10,'m','LineWidth',2)
            plot(t,mean(DNGoodStorm),'k','LineWidth',1)
            %plot(tc,DDensityStorm1Day,'k.','MarkerSize',30)
            title(sprintf('GOES-%d; %s-%s',sat,datestr(Io),datestr(If)))
            xlabel('Time since onset [hrs]')
            legend('Density [amu/cm^3]','-Dst [nT]',...
                    '-Dst [nT] (OMNI)','F10.7/10',...
                    '# values','Location','NorthWest')
            fname = sprintf('Dst_Events_GOES%d_%s_%s',...
                    sat,datestr(Io,29),datestr(If,29));
            title(sprintf('GOES-%d; %s-%s; %d Dst Events',...
                    sat,datestr(Io),datestr(If),length(DDstStorm)));
            %set(gca,'XTick',[-24:2:24])
            plotcmds(fname,1)
    end

    I = find(DDensityMedian(1:end-1) < 30 & DDensityMedian(2:end) > 30);
    clear *Storm*

    if (length(I) > 0)
        k = 1;
        for i = 1:length(I)
            a = I(i)-24*Nd;
            b = I(i)+24*Nd;
            if (a < 1),continue,end
            if (b > length(DDensityMedian)),break,end
            DDensityStorm(k,:) = DDensityMedian(a:b);
            DF107Storm(k,:)  = DF107Median(a:b);
            DDstStorm(k,:)  = DDstMedian(a:b);
            DNGoodStorm(k,:) = DNGood(a:b);
            k = k+1;
        end

        t = [-24*Nd:24*Nd];
        figure(2);clf;hold on;grid on;
            plot(t,nanmean(DDensityStorm),'b','LineWidth',3)
            plot(t,nanmean(-DDstStorm),'g','LineWidth',3)
            plot(t,nanmean(DF107Storm)/10,'m','LineWidth',2)
            plot(t,nanmean(DNGoodStorm),'k','LineWidth',1)
            title(sprintf('GOES-%d; %s-%s',sat,datestr(Io),datestr(If)))
            xlabel('Time since onset [hrs]')
            legend('Density [amu/cm^3]','-Dst [nT]',...
                    'F10.7/10','# values','Location','NorthWest')
            fname = sprintf('Density_Events_GOES%d_%s_%s',...
                    sat,datestr(Io,29),datestr(If,29));
            title(sprintf('GOES-%d; %s-%s; %d \\rho_{eq} Events',...
                    sat,datestr(Io),datestr(If),length(I)));
            %set(gca,'XTick',[-24:2:24])
            plotcmds(fname,1)
    end

end
