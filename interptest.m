function interptest(t,x,tnew,dt)

%{
tic
for i=1:length(tnew)
    datanew(i)=mean(x(find(t>(tnew(i)-dt) & t<(tnew(i)+dt)))); 
end
toc
datanew;
%}
tic
t=t-t(1);
tnew=tnew-tnew(1);
t=floor(t.*10000)./10000;
timestep = min(diff(t)); %// Minimum time-stepsize for t
%t_all=t(1):(t(2)-t(1)):t(end); %But this has rounding errors
t_all=linspace(t(1),t(end),437146);
%{
for i=2:length(t_all)
    t_all(i)=addtodate(t_all(i-1),10,'minute');
end
%}
%t_all = min(t):timestep:max(t); %// create all the timesteps
t_all=floor(t_all.*10000)./10000;
[b1,b2] = ismember(t,t_all);

ind = bsxfun(@plus,[tnew(1)-dt:1:tnew(1)+dt]',[0:numel(tnew)-1]*dt);
[v1,v2] = ismember(ind,t_all(b2));
vind = v2~=0;
v2(v2==0) = NaN;
v2(vind) = x(v2(vind));
out = nanmean(v2);

toc

stuff=1;

%{
tic
ttemp=t(1):(t(2)-t(1)):t(end); %But this has rounding errors
for i=2:length(ttemp)
    ttemp(i)=addtodate(ttemp(i-1),10,'minute');
end

datatemp=ones(1,length(ttemp))*NaN;
mymap=containers.Map(t,x);
for i = 1:length(ttemp)
   if(mymap.isKey(ttemp(i)))
       datatemp(i)=mymap(ttemp(i));
   end
end

datanew2=block_mean_nonflag(datatemp,6);
tnew=block_mean_nonflag(ttemp,6);


toc
%}

stuff=1;


%{
tic
tnew_lb = tnew-dt; %// lower bound
tnew_ub = tnew+dt; %// upper bound
[r,c] = find(bsxfun(@gt,t',tnew_lb) & bsxfun(@lt,t',tnew_ub));
datanew = accumarray(c,x(r),[], @mean);
toc
datanew;
%}