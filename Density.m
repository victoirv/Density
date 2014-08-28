close all;clear all;
%load('../OMNI_OMNI2_merged')
%{
getyears=Year>=2000;
getyears=getyears+(Year<=2011);
getyears=(getyears==2);
VBS=VBS(getyears);
DST=Dst_index(getyears);
ION=Ion_density(getyears);
HOUR=Hour(getyears);
KP=Kp_index(getyears);
%}
FILLED=dlmread('WGhourFS_72_13.txt',',',1,0);
headers=textread('WGhourFS_72_13.txt','%s',28,'delimiter',',');
VBS=1/2*FILLED(:,6).*(abs(FILLED(:,5))-FILLED(:,5));
OMNIDensity=FILLED(:,7);
OMNITime=datenum(FILLED(:,1),1,0)+FILLED(:,2)+FILLED(:,3)/24;
%VBS=1/2*Plasma_bulk_speed.*(abs(Bz_GSM)-Bz_GSM);
%OMNITime=Time;

if(exist('DentonDensityAndTime.mat','file'))
    load('DentonDensityAndTime')
else
    IN=dlmread('massdensity.txt');
    IN(IN(:,1)~=6,:)=[]; %Only use satellite 7 (6 should also work)
    %IN=sortrows(IN,2);
    IN(IN==9999)=NaN;

    [DentonTime,uniquerows]=unique(datenum(IN(:,2),1,0) + IN(:,3)+IN(:,4)/24+IN(:,5)/(24*60),'rows');
    MassDensity=IN(uniquerows,88);

    save('DentonDensityAndTime','DentonTime','MassDensity');
end

if(~exist('figures','dir'))
   mkdir('figures'); 
end

gettime=OMNITime>min(DentonTime);
gettime=gettime+(OMNITime<max(DentonTime));
gettime=(gettime==2);

OMNITime=OMNITime(gettime);
VBS=VBS(gettime);
OMNIDensity=OMNIDensity(gettime);
FILLED=FILLED(gettime,:);

MassDensitySpline=interp1(DentonTime,MassDensity,OMNITime,'linear');

N=50;
Na=0;
Nb=1;
lag=0;
advance=0;

x=MassDensitySpline;

if(~exist(sprintf('figures/OMNI_%s.png',headers{1}),'file'))
    
    for i=1:length(headers)  
        close all;
    figure;plot(OMNITime,FILLED(:,i))
    datetick
    ylabel(headers{i})
    xlabel('Time')
    title(sprintf('OMNI %s',headers{i}))
    filestring=sprintf('figures/OMNI_%s.png',headers{i});
    print('-dpng',filestring)
    end

end

corrs_1=1:(length(headers)+1);
Nb=1;
for i=1:length(headers)  
    f=FILLED(:,i);
    [ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    corrs_1(i)=corr;
    %fprintf('%s - %2.5f\n',headers{i},corr);
    
    close all;
    figure; plot(OMNITime,x); hold on; plot(OMNITime,xnew,'r'); legend('Data','Prediction','Location','SouthEast')
ylabel('Density')
xlabel('Time')
datetick
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.5f',headers{i},Na,Nb,corr))

filestring=sprintf('figures/density_%s_%d_%d.png',headers{i},Na,Nb);
print('-dpng',filestring)

figure;plot(-advance:Nb-advance-1,flipud(cb))
grid
ylabel('Impulse Coefficient')
xlabel('Time Lags')
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.5f',headers{i},Na,Nb,corr))
filestring=sprintf('figures/density_%s_%d_%d_Cb.png',headers{i},Na,Nb);
print('-dpng',filestring)
    
end

%Add VBS_1
x=MassDensitySpline;
f=VBS;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
corrs_1(end)=corr;
    close all;
    figure; plot(OMNITime,x); hold on; plot(OMNITime,xnew,'r'); legend('Data','Prediction','Location','SouthEast')
ylabel('Density')
xlabel('Time')
datetick
title(sprintf('Denton Density from OMNI VBS: Nx:%d Nf:%d, corr: %2.2f',Na,Nb,corr))
filestring=sprintf('figures/density_VBS_%d_%d.png',Na,Nb);
print('-dpng',filestring)

figure;plot(-advance:Nb-advance-1,flipud(cb))
grid
ylabel('Impulse Coefficient')
xlabel('Time Lags')
title(sprintf('Denton Density from OMNI VBS: Nx:%d Nf:%d, corr: %2.2f',Na,Nb,corr))
filestring=sprintf('figures/density_VBS_%d_%d_Cb.png',Na,Nb);
print('-dpng',filestring)



%Open the table for writing
table=fopen('table.txt','w');
Nb=120;
fprintf(table,'Variable \t corr(1) \t corr(120)\n');
for i=1:length(headers)
    f=FILLED(:,i);
    [ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    fprintf(table,'%s \t %2.5f \t %2.5f\n',headers{i},corrs_1(i),corr);
    
    close all;
    figure; plot(OMNITime,x); hold on; plot(OMNITime,xnew,'r'); legend('Data','Prediction','Location','SouthEast')
ylabel('Density')
xlabel('Time')
datetick
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.5f',headers{i},Na,Nb,corr))

filestring=sprintf('figures/density_%s_%d_%d.png',headers{i},Na,Nb);
print('-dpng',filestring)

figure;plot(-advance:Nb-advance-1,flipud(cb))
grid
ylabel('Impulse Coefficient')
xlabel('Time Lags')
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.5f',headers{i},Na,Nb,corr))
filestring=sprintf('figures/density_%s_%d_%d_Cb.png',headers{i},Na,Nb);
print('-dpng',filestring)
    
end

%Add VBS 120
x=MassDensitySpline;
f=VBS;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
fprintf(table,'VBS \t %2.5f \t %2.5f\n',corrs_1(end),corr);
fclose(table);

    close all;
    figure; plot(OMNITime,x); hold on; plot(OMNITime,xnew,'r'); legend('Data','Prediction','Location','SouthEast')
ylabel('Density')
xlabel('Time')
datetick
title(sprintf('Denton Density from OMNI VBS: Nx:%d Nf:%d, corr: %2.2f',Na,Nb,corr))

filestring=sprintf('figures/density_VBS_%d_%d.png',Na,Nb);
print('-dpng',filestring)

figure;plot(-advance:Nb-advance-1,flipud(cb))
grid
ylabel('Impulse Coefficient')
xlabel('Time Lags')
title(sprintf('Denton Density from OMNI VBS: Nx:%d Nf:%d, corr: %2.2f',Na,Nb,corr))
filestring=sprintf('figures/density_VBS_%d_%d_Cb.png',Na,Nb);
print('-dpng',filestring)

%{
plot(DentonTime,MassDensity)
hold on; plot(OMNITime,OMNIDensity,'r')
title('Denton Density vs OMNI corrected Density')
datetick
print -dpng densitycomp.png

MassDensitySpline=spline(DentonTime,MassDensity,OMNITime);

plot(DentonTime,MassDensity)
hold on; plot(OMNITime,MassDensitySpline,'r')
title('Denton Density vs Splined Denton Density')
datetick
print -dpng densitysplinedensity.png

ccoef=corrcoef(MassDensitySpline,OMNIDensity);
fprintf('Correlation between Denton spline and OMNI density: %2.2f',ccoef(1,2))

N=50;
Na=0;
Nb=48;
lag=0;
advance=16;

x=MassDensitySpline;
f=OMNIDensity;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);

figure; plot(x); hold on; plot(xnew,'r'); legend('Data','Prediction','Location','SouthEast')
ylabel('Density response to Density')
xlabel('Time Lags')
title(sprintf('Denton Density from OMNI Density: Nx:%d Nf:%d, corr: %2.2f',Na,Nb,corr))

filestring=sprintf('density_density_%d_%d.png',Na,Nb);
print('-dpng',filestring)

figure;
plot(flipud(ca));

%--------------------------
x=MassDensitySpline;
f=VBS;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);

figure; plot(x); hold on; plot(xnew,'r'); legend('Data','Prediction','Location','SouthEast')
ylabel('Density response to VBs')
xlabel('Time Lags')
title(sprintf('Denton Density from OMNI VBS: Nx:%d Nf:%d, corr: %2.2f',Na,Nb,corr))

filestring=sprintf('density_vbs_%d_%d.png',Na,Nb);
print('-dpng',filestring)


%--------------------------
x=OMNIDensity;
f=VBS;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);

figure; plot(x); hold on; plot(xnew,'r'); legend('Data','Prediction','Location','SouthEast')
ylabel('Density response to VBs')
xlabel('Time Lags')
title(sprintf('OMNI Density from OMNI VBS: Nx:%d Nf:%d, corr: %2.2f',Na,Nb,corr))

filestring=sprintf('density_vbs_%d_%d.png',Na,Nb);
print('-dpng',filestring)

%Do Denton Density for all OMNI variables

%}
