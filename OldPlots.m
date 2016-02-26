%%%%%%%%%%%%%%%%%%%%%%
%Older plots, possibly for thesis, just batch comparing variables and
%making correlation tables



VBS=1/2*FILLED(:,6).*(abs(FILLED(:,5))-FILLED(:,5));
VBz=FILLED(:,6).*FILLED(:,5);
BS=1/2*(abs(FILLED(:,5))-FILLED(:,5));
OMNIDensity=FILLED(:,7);

%Get re-scaled values for the sake of plotting overlays
OverlayFilled=FILLED;
OHr=FILLED(:,3);
DHr=hour(DentonTime);
for i=1:length(headers)
    OverlayFilled(:,i)=((250-150).*(OverlayFilled(:,i)-min(OverlayFilled(:,i))))./(max(OverlayFilled(:,i))-min(OverlayFilled(:,i))) + 150;
end
OverlayBS=((250-150).*(BS-min(BS)))./(max(BS)-min(BS)) + 150;
OverlayVBS=((250-150).*(VBS-min(VBS)))./(max(VBS)-min(VBS)) + 150;
OverlayVBz=((250-150).*(VBz-min(VBz)))./(max(VBz)-min(VBz)) + 150;
OverlayF107=((250-150).*(F107-min(F107)))./(max(F107)-min(F107)) + 150;
OverlayDBz=((250-150).*(DBz-min(DBz)))./(max(DBz)-min(DBz)) + 150;
OverlayDBS=((250-150).*(DBS-min(DBS)))./(max(DBS)-min(DBS)) + 150;


if(MakePlots) %Stack plots
    for i=1:length(headers)
        subplot('position',subplotstack(5,mod(i,4)+1));plot(AVs(:,i),'+-');
        ylabel(headers{i})
        if(mod(i,4)==0 || i==length(headers))
            subplot('position',subplotstack(5,5)); [AX,H1,H2]=plotyy(1:length(AVMDs),AVMDs(1,:),1:length(AVMDs),AVnnans,'plot','bar');
            set(AX(2),'Xlim',[1 length(AVMDs)]); set(AX(1),'Xlim',[1 length(AVMDs)]); set(get(H2,'child'),'facea',.3); set(H1,'marker','+','color','red'); set(AX(1),'YColor','r');
            %drawnow; hold on; xlim manual; h=bar(AVnnans); ch = get(h,'child'); set(ch,'facea',.3)
            ylabel(AX(1),'\rho_{eq}'); ylabel(AX(2),'# of values');
            %xlabel('Time before storm onset (hr)')
            print('-dpng',sprintf('figures/stormavs-%d',i));
            close all
        end
    end
end

if(MakePlots)
    NanPH=zeros(1,24);
    for i=0:23
        NanPH(i+1)=sum(isnan(MassDensitySpline(hour(FILLEDTime)==i)));
    end
    plot(0:23,NanPH);
    ylabel('Nans per hour')
    xlabel('Hour (UT)')
    print -dpng figures/nanph.png
end

if(MakePlots) %Make plot showing that most of F10.7 and Mass Density correlation is from long term effects
    [~, cb, ~,xnew,corr] = IR(log(MassDensity),log(F107),0,12,0,0);
    h=figure('Visible',visible);
    plot(DentonTime,log(MassDensity));
    hold on;
    plot(DentonTime,xnew,'g')
    xlabel('Time')
    ylabel('log(Mass Density)')
    legend('Data','Predicted')
    title(sprintf('Predicting Mass Density from F10.7 - cc(12)=%2.2f',corr))
    print -dpng figures/longtermf107.png
end

if(MakePlots)
    plot(OverlayBS)
    hold on; plot(OverlayVBz,'r')
    plot(MassDensitySpline,'g')
    plot(OverlayFilled(:,15),'k')
    legend('BS','VBz','MassDensity','DST','F10.7','Location','SouthWest')
    timestr=datestr(startt)
    title(timestr)
    print -dpng figures/onestorm-4.png
end


Na=0;
lag=0;
advance=0;
slice=0;
slicemin=10;
slicemax=24;

%x=detrend(MassDensitySpline);
x=MassDensitySpline;

%{
%Figure out what's going on with the year
[ca, cb, cc,xnew,corr] = IR(x,FILLED(:,1),Na,Nb,lag,advance);
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

%Calculate correlation coefficients
cc=corrcoef(F107,DVsw);
fprintf('F10.7 Vsw: %2.2f \n',cc(2))
cc=corrcoef(F107,DBz);
fprintf('F10.7 Bz: %2.2f \n',cc(2))
cc=corrcoef(x,FILLED(:,15),'rows','complete');
fprintf('DST MassDensity: %2.2f \n',cc(2))

if(MakePlots)
    for i=0:23
        avf107(i+1)=mean(F107(DHr==i));
        avdens(i+1)=mean(MassDensity(DHr==i));
        avvsw(i+1)=mean(DVsw(DHr==i));
        avbz(i+1)=mean(DBz(DHr==i));
    end
    avf107=(avf107-min(avf107))/(max(avf107)-min(avf107));
    avdens=(avdens-min(avdens))/(max(avdens)-min(avdens));
    avvsw=(avvsw-min(avvsw))/(max(avvsw)-min(avvsw));
    avbz=(avbz-min(avbz))/(max(avbz)-min(avbz));
    h=figure('Visible',visible);
    plot(0:23,avf107,'r+-')
    hold on
    plot(0:23,avdens,'b+-')
    plot(0:23,avvsw,'k+-')
    plot(0:23,avbz,'m+-')
    legend('F107','Density','Vsw','Bz')
    xlabel('UT Hour')
    ylabel('Normalized Average')
    print '-dpng' 'figures/avf107.png'
end

%Add extra variables
headers=[headers;{'VBS';'BS'; 'VBz'}];
FILLED=[FILLED, VBS, BS, VBz];
OverlayFilled=[OverlayFilled, OverlayVBS, OverlayBS, OverlayVBz];

dheaders={'DBS';'F107';'lnF107';'DBz'};
DFILLED=[DBS, F107, log(F107), DBz];
DOverlayFilled=[OverlayDBS, OverlayF107, log(OverlayF107), OverlayDBz];

if(MakePlots && ~exist(sprintf('figures/OMNI_%s.png',headers{1}),'file'))
    for i=1:length(headers)
        close all;
        h=figure('Visible',visible);plot(FILLEDTime,FILLED(:,i))
        datetick
        ylabel(headers{i})
        xlabel('Time')
        title(sprintf('OMNI %s',headers{i}))
        filestring=sprintf('figures/OMNI_%s.png',headers{i});
        print('-dpng',filestring)
    end
end
if(MakePlots && ~exist(sprintf('figures/Denton_%s.png',dheaders{1}),'file'))
    for i=1:length(dheaders)
        close all;
        h=figure('Visible',visible);plot(DentonTime,DFILLED(:,i))
        datetick
        ylabel(dheaders{i})
        xlabel('Time')
        title(sprintf('Denton %s',dheaders{i}))
        filestring=sprintf('figures/Denton_%s.png',dheaders{i});
        print('-dpng',filestring)
    end
end
if(MakePlots && ~exist(sprintf('figures/Density_vs_%s.png',dheaders{1}),'file'))
    for i=1:length(headers)
        close all;
        h=figure('Visible',visible);plot(MassDensitySpline,FILLED(:,i),'+')
        ylabel(headers{i})
        xlabel('Mass Density')
        cc=corrcoef(MassDensitySpline,FILLED(:,i),'rows','complete');
        title(sprintf('%s - cc:%2.2f',headers{i},cc(2)))
        filestring=sprintf('figures/Density_vs_%s.png',headers{i});
        print('-dpng',filestring)
    end
end

%x=detrend(MassDensitySpline);

%Open the table for writing
if(slice)
    table=fopen('tableslice.txt','w');
else
    table=fopen('table.txt','w');
end
fprintf(table,'<pre>\n');
Nb1=1; %Just in case?
Nb=12;
x=log(MassDensitySpline);
fprintf(table,'Input \t CC1 \t CC12 \t PE(1) \t PE(12)\n');
BigTable={};
for i=1:length(headers)
    f=FILLED(:,i);
    if(slice)
        [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,OHr,slicemin,slicemax);
        [~,cb,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,OHr,slicemin,slicemax);
    else
        [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
        [~, cb, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
    end
    pe1=pe_nonflag(x,xnew1);
    if(MakePlots), densityplot(FILLEDTime,[x,xnew1,OverlayFilled(:,i)],headers{i},Na,Nb1,corr1,pe1,visible), end
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{headers{i},corr1,corr,pe1,pe}];
    if(MakePlots), densityplot(FILLEDTime,[x,xnew,OverlayFilled(:,i)],headers{i},Na,Nb,corr,pe,visible), end
    if(MakePlots), densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe,visible), end
end

x=log(MassDensity);
for i=1:length(dheaders)
    f=DFILLED(:,i);
    if(slice)
        [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,DHr,slicemin,slicemax);
        [~,cb,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,DHr,slicemin,slicemax);
    else
        [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
        [~, cb, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
    end
    pe1=pe_nonflag(x,xnew1);
    if(MakePlots), densityplot(DentonTime,[x,xnew1,DOverlayFilled(:,i)],dheaders{i},Na,Nb1,corr1,pe1,visible), end
    pe=pe_nonflag(x,xnew);
    BigTable=[BigTable;{dheaders{i},corr1,corr,pe1,pe}];
    if(MakePlots), densityplot(DentonTime,[x,xnew,DOverlayFilled(:,i)],headers{i},Na,Nb,corr,pe,visible), end
    if(MakePlots), densitycoefplot(-advance:Nb-advance-1,flipud(cb),headers{i},Na,Nb,corr,pe,visible), end
end

%Now for doubles
x=log(MassDensitySpline);
f=[FILLED(:,3) FILLED(:,5)]; %Hr and Bz
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,OHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,OHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'Hr+Bz',corr1,corr,pe1,pe}];

f=[FILLED(:,5) FILLED(:,6)]; %Bz and V
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,OHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,OHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'Bz+V',corr1,corr,pe1,pe}];

x=log(MassDensity);
f=[F107 DBz]; %F107 and Bz
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,DHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,DHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'F107+Bz',corr1,corr,pe1,pe}];
f=[F107 DVsw]; %Bz and V
if(slice)
    [~,~,~,xnew1,corr1] = IRslice(x,f,Na,Nb1,lag,advance,DHr,slicemin,slicemax);
    [~,~,~,xnew,corr] = IRslice(x,f,Na,Nb,lag,advance,DHr,slicemin,slicemax);
else
    [~,~,~,xnew1,corr1] = IR(x,f,Na,Nb1,lag,advance);
    [~, ~, ~,xnew,corr] = IR(x,f,Na,Nb,lag,advance);
end
pe1=pe_nonflag(x,xnew1);
pe=pe_nonflag(x,xnew);
BigTable=[BigTable;{'F107+V',corr1,corr,pe1,pe}];


BigTable=sortrows(BigTable,2);
for i=1:size(BigTable,1)
    if(BigTable{i,5}>-10)
        fprintf(table,'%s\t %2.2f \t %2.2f \t %2.2f \t %2.2f\n',BigTable{i,1},BigTable{i,2},BigTable{i,3},BigTable{i,4},BigTable{i,5});
    else
        fprintf(table,'%s\t %2.2f \t %2.2f \t %2.2f \t %2.0e\n',BigTable{i,1},BigTable{i,2},BigTable{i,3},BigTable{i,4},BigTable{i,5});
    end
end
fprintf(table,'\n</pre>');

fclose(table);
system('cat README.txt table.txt > README.md');


end

function densityplot(x,ys,string,Na,Nb,corr,eff,visible)
close all;
test=size(ys);
if test(1)>test(2)
    ys=ys';
end
h=figure('Visible',visible);plot(x,ys)
legend('Data','Prediction',string,'Location','NorthEast')
ylabel('Density')
xlabel('Time')
datetick
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
filestring=sprintf('figures/density_%s_%d_%d.png',string,Na,Nb);
print(h,'-dpng',filestring)
end

function densitycoefplot(x,ys,string,Na,Nb,corr,eff,visible)
close all;
h=figure('Visible',visible);
plot(x,ys)
ylabel('Impulse Coefficient')
xlabel('Time Lags')
grid
title(sprintf('Denton Density from OMNI %s: Nx:%d Nf:%d, corr: %2.3f, eff: %2.3f',string,Na,Nb,corr,eff))
filestring=sprintf('figures/density_%s_%d_%d_Cb.png',string,Na,Nb);
print(h,'-dpng',filestring)
end
