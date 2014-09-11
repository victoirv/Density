function Density

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
BS=abs(FILLED(:,5))-FILLED(:,5);
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
OverlayFilled=FILLED;
%Rescale to overlay on plots of density (where there's space from 150-250)
for i=1:length(headers)
    OverlayFilled(:,i)=((250-150).*(OverlayFilled(:,i)-min(OverlayFilled(:,i))))./(max(OverlayFilled(:,i))-min(OverlayFilled(:,i))) + 150;
end
    
MassDensitySpline=interp1(DentonTime,MassDensity,OMNITime,'linear');



%Compare the two densities
figure; plot(OMNITime,MassDensitySpline);hold on; plot(OMNITime,OMNIDensity,'r')
legend('Denton','OMNI','Location','NorthEast')
title('OMNI Density vs Denton Density')
ylabel('Density')
xlabel('Time')
datetick
print -dpng figures/densitycomp.png

Na=0;
Nb=1;
lag=0;
advance=0;

x=detrend(MassDensitySpline);

%{ 
%Figure out what's going on with the year
[ca, cb, cc,xnew,corr] = IRboot(x,FILLED(:,1),Na,Nb,lag,advance);
figure; plot(xnew,MassDensitySpline,'+')

xm=mean(x);
ym=mean(FILLED(:,1));
sx=std(x);
sy=std(FILLED(:,1));
zx=(x-xm)./sx;
zy=(FILLED(:,1)-ym)./sy;
c=sum((zx.*zy))/(length(zx)-1);
title(sprintf('Denton density vs Year predicted density, my c: %2.2f - corrcoef c: %2.2f',c,corr));
print -dpng figures/yearmystery.png

lx=length(x);
lx5=floor(length(x)/5);
test=ones(1,lx);
test(lx5:2*lx5)=2;
test(2*lx5:3*lx5)=3;
test(3*lx5:4*lx5)=4;
test(4*lx5:5*lx5)=5;

test2=1:lx;

corrcoef(test,x)
%}

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
pes_1=1:(length(headers)+1);
Nb=1;
for i=1:length(headers)
    f=FILLED(:,i);
    [ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    corrs_1(i)=corr;  
    pe=pe_nonflag(x,xnew);
    pes_1(i)=pe;
    densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe)
    
end

%Add VBS_1
x=detrend(MassDensitySpline);
f=VBS;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
corrs_1(end)=corr;
pes_1(end)=pe_nonflag(x,xnew);

densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],'VBS',Na,Nb,corr,pe_nonflag(x,xnew))

%Open the table for writing
table=fopen('table.txt','w');
Nb=120;
fprintf(table,'Variable \t corr(1) \t corr(120) \t eff(1) \t eff(120)\n');
for i=1:length(headers)
    f=FILLED(:,i);
    [ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
    pe=pe_nonflag(x,xnew);
    fprintf(table,'%s \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',headers{i},corrs_1(i),corr,pes_1(i),pe);
    
    densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe)
    densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe)
end

%Add VBS 120
f=VBS;

[ca, cb, cc,xnew,corr] = IRboot(x,f,Na,Nb,lag,advance);
pe=pe_nonflag(x,xnew);
fprintf(table,'VBS \t %2.5f \t %2.5f \t %2.5f \t %2.5f\n',corrs_1(end),corr,pes_1(end),pe);
fclose(table);
system('cat README.txt table.txt > README.md');

densityplot(OMNITime,[x,xnew,OverlayFilled(:,i)],'VBS',Na,Nb,corr,pe)
densitycoefplot(-advance:Nb-advance-1,flipud(cb),'VBS',Na,Nb,corr,pe)
end

function densityplot(x,ys,string,Na,Nb,corr,eff)
close all;
test=size(ys);
if test(1)>test(2) 
   ys=ys'; 
end
figure;plot(x,ys)
    legend('Data','Prediction',string,'Location','NorthEast')
    ylabel('Density')
    xlabel('Time')
    datetick
    title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
    filestring=sprintf('figures/density_%s_%d_%d.png',string,Na,Nb);
print('-dpng',filestring)
end

function densitycoefplot(x,ys,string,Na,Nb,corr,eff)
close all;
figure;plot(x,ys)
ylabel('Impulse Coefficient')
xlabel('Time Lags')
grid    
    title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
    filestring=sprintf('figures/density_%s_%d_%d_Cb.png',string,Na,Nb);
    print('-dpng',filestring)
end