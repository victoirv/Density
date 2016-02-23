function NNAnalysis(AVMat,AVMDMat)

loops=40;
disp('NN - Bz')
x=nanmedian(AVMat(:,20:24,30),2);
z=nanmedian(AVMDMat(:,25:29),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'F10.7','rho_eq'},satnum);

x=nanmedian(AVMat(:,20:24,5),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'Bz','rho_eq'},satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[5 30]),2),z,1,loops,{'Bz','F10.7','rho_eq'},satnum);


%Same thing but Vsw
disp('NN - Vsw')
x=nanmedian(AVMat(:,20:24,6),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'Vsw','rho_eq'},satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[6 30]),2),z,1,loops,{'Vsw','F10.7','rho_eq'},satnum);

%Same thing but MLT
disp('NN - MLT')
xtest=linspace(0,24);
x=nanmedian(AVMat(:,20:24,29),2);
nnbehavior(x,z,xtest,1,loops,{'MLT','rho_eq'},satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[29 30]),2),z,1,loops,{'MLT','F10.7','rho_eq'},satnum);

%DST
disp('NN - Dst')
x=nanmedian(AVMat(:,20:24,15),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'Dst','rho_eq'},satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[15 30]),2),z,1,loops,{'Dst','F10.7','rho_eq'},satnum);

%Doy
disp('NN - Doy')
x=nanmedian(AVMat(:,20:24,2),2);
xtest=linspace(1,365);
nnbehavior(x,z,xtest,1,loops,{'Doy','rho_eq'},satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[2 30]),2),z,1,loops,{'Doy','F10.7','rho_eq'},satnum);

%Rhosw
disp('NN - Rhosw')
x=nanmedian(AVMat(:,20:24,31),2);
xtest=linspace(min(x),max(x));
nnbehavior(x,z,xtest,1,loops,{'Rhosw','rho_eq'},satnum);

nnbehavior2(nanmedian(AVMat(:,20:24,[31 30]),2),z,1,loops,{'Rhosw','F10.7','rho_eq'},satnum);

%Print correlations as table of linear vs neural net coefficients
fprintf('Figures done. Doing table\n');
PermNames={'doy','Bz','Vsw','Dst','MLT','F107','Rhosw'};
PermCols=[2 5 6 15 29 30 31];

table=fopen(sprintf('tables/NNtable-GOES%d.txt',satnum),'w');
fprintf(table,'<pre>\n');
fprintf(table,'Vars \t \t  CCtr  NNtr  CCt   NNt   CCv   NNv\n');
Perms=combnk(1:7,1);
for i=1:length(Perms)
    CCMs(:,:)=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
    fprintf(table,'%s      \t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:));
end
Perms=combnk(1:7,2);
for i=1:length(Perms)
    CCMs(:,:)=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
    fprintf(table,'%s   \t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:));
end
Perms=combnk(1:7,3);
for i=1:length(Perms)
    CCMs(:,:)=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(i,:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
    fprintf(table,'%s   \t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f \n',strjoin(PermNames(Perms(i,:)),'+'),CCMs(:));
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
CCMs(:,:)=nntest(nanmedian(AVMat(:,20:24,PermCols(Perms(:))),2),nanmedian(AVMDMat(:,25:29),2),1,loops,1);
fprintf(table,'%s\t- %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f %+2.2f \n',strjoin(PermNames(Perms(:)),'+'),CCMs(:));

fclose(table);