function CCAnalysis(AVMat,AVMDMat,satnum)
loops=100;
trainpercent=0.50;

table=fopen(sprintf('tables/CCtable-GOES%d.txt',satnum),'w');
fprintf(table,'<pre>\n');
fprintf(table,'Vars \t \t  CCtr  CCt +-CCtr +-CCt\n');

PermNames={'DoY','B_z','V_{sw}','D_{st}','MLT','F_{10.7}','\rho_{sw}'};
PermCols=[2 5 6 15 29 30 31];

Perms=combnk(1:7,1);
for i=1:length(Perms)
    x=nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2);
    y=nanmedian(AVMDMat(:,25:29),2);
    Z=[x y];
    for col=1:min(size(Z))
        Z(isnan(Z(:,col)),:)=[];
    end
    x=Z(:,1:end-1);
    y=Z(:,end);
    
    for j=1:loops
        tri=randsample(1:length(y),floor(trainpercent*length(y)));
        ti=setdiff(1:length(y),tri);
        
        Ztr=[x(tri,:) ones(length(tri),1) y(tri)];
        Zt=[x(ti,:) ones(length(ti),1) y(ti)];
        coef=Ztr(:,1:end-1)\Ztr(:,end);
        cctrm=corrcoef(Ztr(:,1:end-1)*coef,Ztr(:,end));
        cctm=corrcoef(Zt(:,1:end-1)*coef,Zt(:,end));
        cct(j,1)=cctrm(1,2);
        cct(j,2)=cctm(1,2);
    end
    CCMs=nanmedian(cct);
    CCMsds=nanstd(cct);
    fprintf(table,'%s      \t- %+2.2f %+2.2f %2.2f %2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:),CCMsds(:));
end

Perms=combnk(1:7,2);
for i=1:length(Perms)
    x=nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2);
    x=reshape(x,length(x),[]);
    y=nanmedian(AVMDMat(:,25:29),2);
    Z=[x y];
    for col=1:min(size(Z))
        Z(isnan(Z(:,col)),:)=[];
    end
    x=Z(:,1:end-1);
    y=Z(:,end);
    
    for j=1:loops
        tri=randsample(1:length(y),floor(trainpercent*length(y)));
        ti=setdiff(1:length(y),tri);
        
        Ztr=[x(tri,:) ones(length(tri),1) y(tri)];
        Zt=[x(ti,:) ones(length(ti),1) y(ti)];
        coef=Ztr(:,1:end-1)\Ztr(:,end);
        cctrm=corrcoef(Ztr(:,1:end-1)*coef,Ztr(:,end));
        cctm=corrcoef(Zt(:,1:end-1)*coef,Zt(:,end));
        cct(j,1)=cctrm(1,2);
        cct(j,2)=cctm(1,2);
    end
    CCMs=nanmedian(cct);
    CCMsds=nanstd(cct);
    fprintf(table,'%s   \t- %+2.2f %+2.2f %2.2f %2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:),CCMsds(:));
end

Perms=combnk(1:7,7);
x=nanmedian(AVMat(:,20:24,PermCols(Perms(:))),2);
x=reshape(x,length(x),[]);
y=nanmedian(AVMDMat(:,25:29),2);
Z=[x y];
for col=1:min(size(Z))
    Z(isnan(Z(:,col)),:)=[];
end
x=Z(:,1:end-1);
y=Z(:,end);

for j=1:loops
    tri=randsample(1:length(y),floor(trainpercent*length(y)));
    ti=setdiff(1:length(y),tri);
    
    Ztr=[x(tri,:) ones(length(tri),1) y(tri)];
    Zt=[x(ti,:) ones(length(ti),1) y(ti)];
    coef=Ztr(:,1:end-1)\Ztr(:,end);
    cctrm=corrcoef(Ztr(:,1:end-1)*coef,Ztr(:,end));
    cctm=corrcoef(Zt(:,1:end-1)*coef,Zt(:,end));
    cct(j,1)=cctrm(1,2);
    cct(j,2)=cctm(1,2);
end
CCMs=nanmedian(cct);
CCMsds=nanstd(cct);
fprintf(table,'All\t- %+2.2f %+2.2f %2.2f %2.2f \n',CCMs(:),CCMsds(:));

fclose(table);