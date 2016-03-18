function NNAnalysis(AVMat,AVMDMat,satnum)

loops=40;
AXmap=containers.Map({'B_z','V_{sw}','F_{10.7}','\rho_{eq}','DoY','\rho_{sw}','D_{st}','MLT'},{[-10 10],[200 800],[50 350],[0 70],[0 356],[0 30],[-100 30],[0 24]});


disp('NN - B_z')
x=nanmedian(AVMat(:,20:24,30),2);
z=nanmedian(AVMDMat(:,25:29),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'F_{10.7}','\rho_{eq}'},AXmap,satnum);

x=nanmedian(AVMat(:,20:24,5),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'B_z','\rho_{eq}'},AXmap,satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[5 30]),2),z,1,loops,{'B_z','F_{10.7}','\rho_{eq}'},AXmap,satnum);


%Same thing but V_{sw}
disp('NN - V_{sw}')
x=nanmedian(AVMat(:,20:24,6),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'V_{sw}','\rho_{eq}'},AXmap,satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[6 30]),2),z,1,loops,{'V_{sw}','F_{10.7}','\rho_{eq}'},AXmap,satnum);

%Same thing but MLT
disp('NN - MLT')
xtest=linspace(0,24);
x=nanmedian(AVMat(:,20:24,29),2);
nnbehavior(x,z,xtest,1,loops,{'MLT','\rho_{eq}'},AXmap,satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[29 30]),2),z,1,loops,{'MLT','F_{10.7}','\rho_{eq}'},AXmap,satnum);

%D_{st}
disp('NN - D_{st}')
x=nanmedian(AVMat(:,20:24,15),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'D_{st}','\rho_{eq}'},AXmap,satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[15 30]),2),z,1,loops,{'D_{st}','F_{10.7}','\rho_{eq}'},AXmap,satnum);

%DoY
disp('NN - DoY')
x=nanmedian(AVMat(:,20:24,2),2);
xtest=linspace(1,365);
nnbehavior(x,z,xtest,1,loops,{'DoY','\rho_{eq}'},AXmap,satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[2 30]),2),z,1,loops,{'DoY','F_{10.7}','\rho_{eq}'},AXmap,satnum);

%\rho_{sw}
disp('NN - \rho_{sw}')
x=nanmedian(AVMat(:,20:24,31),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'\rho_{sw}','\rho_{eq}'},AXmap,satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[31 30]),2),z,1,loops,{'\rho_{sw}','F_{10.7}','\rho_{eq}'},AXmap,satnum);

%Print correlations as table of linear vs neural net coefficients
fprintf('Figures done. Doing table\n');
PermNames={'DoY','B_z','V_{sw}','D_{st}','MLT','F_{10.7}','\rho_{sw}'};
PermCols=[2 5 6 15 29 30 31];

table=fopen(sprintf('tables/NNtable-GOES%d.txt',satnum),'w');
fprintf(table,'<pre>\n');
fprintf(table,'Vars \t \t  CCtr  NNtr  CCt   NNt   CCv   NNv  +-CCtr  +-NNtr  +-CCt   +-NNt   +-CCv   +-NNv\n');
Perms=combnk(1:7,1);
for i=1:length(Perms)
    [CCMs(:,:), CCMSDs(:,:)]=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
    fprintf(table,'%s      \t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:),CCMSDs(:));
end
Perms=combnk(1:7,2);
for i=1:length(Perms)
    [CCMs(:,:), CCMSDs(:,:)]=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
    fprintf(table,'%s   \t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:),CCMSDs(:));
end
Perms=combnk(1:7,3);
for i=1:length(Perms)
    [CCMs(:,:), CCMSDs(:,:)]=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
    fprintf(table,'%s   \t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:),CCMSDs(:));
end
%{
    Perms=combnk(1:6,4);
    for i=1:length(Perms)
        CCMs(:,:)=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
        fprintf(table,'%s\t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:));
    end
    Perms=combnk(1:6,5);
    for i=1:length(Perms)
        CCMs(:,:)=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
        fprintf(table,'%s\t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:));
    end
%}
Perms=combnk(1:7,7);
[CCMs(:,:), CCMSDs(:,:)]=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
fprintf(table,'All\t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %2.2f %2.2f %2.2f %2.2f %2.2f %2.2f \n',CCMs(:),CCMSDs(:));

fclose(table);