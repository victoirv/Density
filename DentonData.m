%%%%%%%%%%%%%%%%%%
%Script for reading Denton, Kondrashov, and OMNI data
%Requires "inputs" of satnum, cutyears, and TakahashiCond


%Read filled dataset from Kondrashov(2014)
if(~exist('data/WGhourFS_72_13.txt','file'))
    urlwrite('http://mag.gmu.edu/git-data/victoirv/Density/data/WGhourFS_72_13.txt','data/WGhourFS_72_13.txt');
end

filledname=sprintf('data/FILLED_%d',satnum);
if(TakahashiCond), filledname=sprintf('%s_tak',filledname);end
if(exist(sprintf('%s.mat',filledname),'file'))
    load(filledname)
else
    %Data from supplementary section of http://onlinelibrary.wiley.com/doi/10.1002/2014GL059741/full
    FILLED=dlmread('data/WGhourFS_72_13.txt',',',1,0); 
    headers=textread('data/WGhourFS_72_13.txt','%s',28,'delimiter',',');
    FILLED=FILLED(FILLED(:,1)>1980,:); %Just to get into Denton time
    if(TakahashiCond)
        FILLED=FILLED(FILLED(:,1)>1988,:);
        FILLED=FILLED(FILLED(:,1)<1992,:); %For comparing to Takahashi 2010 Fig 11
    end
    
    OMNIDensity=FILLED(:,7);
    FILLEDTime=datenum(FILLED(:,1),1,0)+FILLED(:,2)+FILLED(:,3)/24;
    
    save(filledname,'FILLED','headers','OMNIDensity','FILLEDTime');
end



%Get Denton data (mostly just for mass density)
if(exist(sprintf('data/DentonDensityAndTime_%d.mat',satnum),'file'))
    load(sprintf('data/DentonDensityAndTime_%d',satnum))
else
    if(~exist('data/massdensity.txt','file'))
        urlwrite('http://mag.gmu.edu/git-data/victoirv/Density/data/massdensity.txt','data/massdensity.txt');
    end
    IN=dlmread('data/massdensity.txt');
    IN(IN(:,1)~=satnum,:)=[]; %Satellites 6 and 7 should work
    IN(IN==9999)=NaN;
    
    [DentonTime,uniquerows]=unique(datenum(IN(:,2),1,0) + IN(:,3)+IN(:,4)/24+IN(:,5)/(24*60),'rows');
    
    %88 is density, 10 is f10.7, 31 is BZ_sw, 32 is Vsw
    %F107=IN(uniquerows,10);
    DBz=IN(uniquerows,31);
    DVsw=IN(uniquerows,32);
    DBS=1/2*IN(uniquerows,32).*(abs(DBz)-DBz);
    AE=IN(uniquerows,12);
    MLT=IN(uniquerows,8);
    MassDensity=IN(uniquerows,88); %88 is rho_eq
    M=IN(uniquerows,88)./IN(uniquerows,68); %65 is nei_n or electron density topside ionosphere
    
    save(sprintf('data/DentonDensityAndTime_%d',satnum),'DentonTime','MassDensity','DBz','MLT','DBS','DVsw','AE','M');
end

%Load F10.7 specifically 
if(exist('data/OMNIdata.mat','file'))
    load('data/OMNIdata.mat')
else
    OMNIin=dlmread('data/omni2_12536.lst');
    OMNIin(OMNIin==999.9)=NaN;
    OMNIin(OMNIin==9999)=NaN;
    OMNITime=datenum(OMNIin(:,1),1,OMNIin(:,2),OMNIin(:,3),0,0);
    F107=OMNIin(:,13); 
    OMNIRho=OMNIin(:,7);
    F107(F107==999.9)=NaN;
    save('data/OMNIdata','OMNITime','F107','OMNIRho');
end

if(~exist('figures','dir'))
    mkdir('figures');
end

%Save interpolated values since it takes a while to interpolate the whole
%dataset. "Spline" is a misnomer at this point but a proper find/replace might
%take a while and introduce bugs
interpname=sprintf('data/InterpedVals_%d',satnum);
if(TakahashiCond), interpname=sprintf('%s_tak',interpname);end
if(exist(sprintf('%s.mat',interpname),'file'))
    load(interpname)
else
    MassDensitySpline=interptest(DentonTime,MassDensity,FILLEDTime,(FILLEDTime(2)-FILLEDTime(1))/2);
    MLTFit=interptest(DentonTime,MLT,FILLEDTime,(FILLEDTime(2)-FILLEDTime(1))/2);
    AEFit=interptest(DentonTime,AE,FILLEDTime,(FILLEDTime(2)-FILLEDTime(1))/2);
    MFit=interptest(DentonTime,M,FILLEDTime,(FILLEDTime(2)-FILLEDTime(1))/2);
    %F107Spline=interptest(DentonTime,F107,FILLEDTime,(FILLEDTime(2)-FILLEDTime(1))/2);
    MassDensitySpline=MassDensitySpline';
    %F107Spline=F107Spline';
    save(interpname,'MassDensitySpline','MLTFit','AEFit','MFit');%,'F107Spline');
end


%Get everything onto the same time grid
gettime=FILLEDTime>min(DentonTime);
gettime=gettime+(FILLEDTime<max(DentonTime));
gettime=gettime+(FILLEDTime<cuttime(2)); %If user wants to cut down time range. 
gettime=gettime+(FILLEDTime>cuttime(1));
gettime=(gettime==4);



%Move everything already created onto desired time grid
FILLEDTime=FILLEDTime(gettime);
OMNIDensity=OMNIDensity(gettime);
FILLED=FILLED(gettime,:);
MassDensitySpline=MassDensitySpline(gettime);
MLTFit=MLTFit(gettime);
AEFit=AEFit(gettime);
MFit=MFit(gettime);

%Start and end year, for titles
sy=str2num(datestr(FILLEDTime(1),10)); ey=str2num(datestr(FILLEDTime(end),10));

%Should only need to match it to FILLEDTime since FILLEDTime is already matched
%to DentonTime at this point
getFtime=OMNITime>=min(FILLEDTime);
getFtime=getFtime+(OMNITime<=max(FILLEDTime));
getFtime=(getFtime==2);
OMNITime=OMNITime(getFtime);
F107=F107(getFtime);
OMNIRho=OMNIRho(getFtime);



%Find storm with enough data to analyze
MassDensityNanSpline=interp1(FILLEDTime(~isnan(MassDensitySpline)),MassDensitySpline(~isnan(MassDensitySpline)),FILLEDTime,'linear');

%FILLED=[FILLED F107Spline]; %When F107 is from Denton, use
%interpolated version
FILLED=[FILLED MLTFit F107 OMNIRho MFit]; %28 cols of F, 29 is mlt, 30 is f107, etc
headers{end+1}='MLT';
headers{end+1}='F_{10.7}';
headers{end+1}='OMNIRho';
headers{end+1}='M';
units={'year','day','hour','nT','nT','cm/s','verify','nPa','','','','','','','nT','nT','nT','nT','nT','nT','nT','','','','','','',''};
units{end+1}='hour';
units{end+1}='s.f.u';
units{end+1}='amu/cm^3';
units{end+1}='amu?';
