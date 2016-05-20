function ModelModels

refpoint=[10 0 0];
run=2;
dataheaders={'x','y','z','bx','by','bz','jx','jy','jz','ux','uy','uz','p','rho'};

if(exist('data/DifferencesData.mat','file'))
    load('data/DifferencesData.mat')
else

    [status, filedata]=system(sprintf('grep -h -e "^%3.6f %3.6f %3.6f" ../Differences/output/Brian_Curtis_042213_%d/data/cuts/*',refpoint(1), refpoint(2), refpoint(3), run)); 
    data=cell2mat(textscan(filedata,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f'));

    inputvars={'Year','Month','Day','Hour','Min','Sec','Msec','Bx[nT]','By[nT]','Bz[nT]','Vx[km/s]','Vy[km/s]','Vz[km/s]','N[cm^(-3)]','T[Kelvin]'};

    inputs=dlmread(sprintf('../Differences/data/Brian_Curtis_042213_%d_IMF.txt',run));

    for i=1:72
        bininputs(i,:)=median(inputs((i-1)*5+1:i*5,:));
    end
    bininputs(73,:)=median(inputs(72*5+1:end,:));
    
    save('data/DifferencesData','data','inputs','bininputs','inputvars');
end

for ctest=8:15
   [~,~,~,~,corr]=IR(data(:,9),bininputs(:,ctest),0,3);
   fprintf('%s: %2.3f\n',inputvars{ctest},corr);
end

[~,~,~,~,corr]=IR(data(:,9),bininputs(:,8:15),0,3);
fprintf('All: %2.3f\n',corr);

plot(data(:,9)); %check for Bz flip

%grep -h -e "-200.000000 0.000000 -47.000000" ../Differences/output/Brian_Curtis_042213_2/data/cuts/*


    test=dlmread('../Differences/output/Brian_Curtis_042213_2/data/cuts/Step_70_Y_eq_0.txt');
    x=unique(test(:,1));
    z=unique(test(:,3));
    d=reshape(test(:,8),length(x),[]);
    surf(x,z,d') %Make sure we're looking at the magnetosphere?
    
    
    %look at columns 8, 10, 11, 12?