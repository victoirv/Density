function [datanew NUsed]=interptest(t,x,tnew,dt)

if nargin<3
    error('datanew=interptest(t,x,tnew,dt)')
end

if nargin == 3
    dt=(tnew(2)-tnew(1))/2; %halfway to the points on either side
end

%{
if(tnew(1)<t(1) || tnew(end)>t(end))
    error('New time must be within existing time')
end
%}

tnew=tnew(:); %Force orientation
t=t(:);
if(size(x,1)~=size(t,1)) %If x has two dimensions and the time dimension isn't aligned with t, flip it (and t should now be a column vector)
    x=x';
end

%{
if(min(size(x))~=1) %If x has two dimensions, recursively interpolate over the non-time dimensions
    xdim=setdiff(size(x),size(t));
    datanew=NaN(length(tnew),xdim);  
    for i=1:xdim
        datanew(:,i)=interptest(t,x(:,i),tnew,dt);
    end
else
    
    while(t(end)>(tnew(end)+dt))
        t(end)=[];
        x(end)=[];
    end

    dt2=floor(dt/(t(2)-t(1)));
    t2=tnew(1)-dt:(t(2)-t(1)):tnew(end)+dt;
    t2(end)=[];  %last point would be first of next bin
    x2=nan(length(t2),1);

    indices=floor((t-t(1))/(t(2)-t(1))+1); %Indices on expanded uniform grid of t (not tnew yet)

    if(tnew(1)>t(1)) %If new time starts later than data provided, cut down x. NEEDS TESTING SINCE RE-WRITE
        offset=find(diff([0 t'>=(tnew(1)-dt) 0])>0)-1;
        x2(indices(offset:end)-indices(offset)+1)=x(offset:end);
    else %New time wants data from before data exists (or tnew(1)==t(1))
        offset=find(diff([0 t(1)<=t2 0])>0)-1; %Find first center point
        x2(indices+offset)=x; %Go from nonuniform time grid to uniform time grid with NaNs
    end

    snip=0;
    while(mod(length(x2),(dt2*2))~=0)
        x2(end+1)=NaN; %Adding nans shouldn't affect final median, just buffer for size
        snip=snip+1;
    end
    datanew=reshape(x2,dt2*2,length(x2)/(dt2*2)); %Reshape and median to smaller uniform time grid
    NUsed=sum(~isnan(datanew))';
    datanew=nanmedian(datanew)';
    
    if(length(datanew)==(length(tnew)+1) && isnan(datanew(end)))
        datanew(end)=[]; %Because sometimes dt2 has a rounding/precision problem and the algorithm buffers up an additional, empty value
    end
    datasave=datanew;

    %{
    if(tnew(1)<t(1))
        prebuffer=tnew(1):(tnew(2)-tnew(1)):t(1);
        prebuffer=prebuffer(1:end-1);
        datanew=[nan(length(prebuffer),1); datanew];
    end
    if(tnew(end)>t(end))
        postbuffer=t(end):(tnew(2)-tnew(1)):tnew(end);
        postbuffer=postbuffer(2:end);
        datanew=[datanew; nan(length(postbuffer,1))];
    end
%}
return

end
%}



xdim=setdiff(size(x),size(t)); %setdiff(size(x),size(t)) finds length of non-time dimension of x
if(isempty(xdim))
    xdim=1;
end
datanew2=zeros(length(tnew),xdim); 

for i=1:length(tnew)
    %datanew(i,:)=nanmedian(x(t>=(tnew(i)-dt) & t<(tnew(i)+dt),:)); %Center
    datanew(i,:)=nanmedian(x(t>=(tnew(i)) & t<(tnew(i)+2*dt),:)); %Forwards time
end
test=0;




%{
tic
%t=t-t(1);
%tnew=tnew-tnew(1);
%t=round(t.*10000)./10000;
timestep = min(diff(t)); %// Minimum time-stepsize for t
%t_all=t(1):(t(2)-t(1)):t(end); %But this has rounding errors
t_all=linspace(t(1),t(end),437146);
%{
for i=2:length(t_all)
    t_all(i)=addtodate(t_all(i-1),10,'minute');
end
%}
t_all = min(t):timestep:t(end); %// create all the timesteps
%t_all=floor(t_all.*10000)./10000;
t = interp1(t_all, t_all, t, 'nearest');
[b1,b2] = ismember(t,t_all);
b2(end)=b2(end-1)+1;
dtnew=min(diff(tnew));
ind = bsxfun(@plus,[linspace(datenum_round_off(tnew(1)-3*dt+dt/2,'minute'),datenum_round_off(tnew(1)+3*dt-dt/2,'minute'),6)]',[0:numel(tnew)-1]*dtnew);
[v1,v2] = ismember(floor(t_all(b2).*10000),floor(ind.*10000));
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
%}