function [starti,endi,duration]=FindStorms(storms,FILLED,cutoffduration,cutconditions,maxwidth)

starti=find(storms>0);
endi=find(storms<0)-1;
duration=endi-starti+1;

%Shift event start points to next local DST minimum 
while(0) 
    ind=FILLED(starti+1,15)<FILLED(starti,15);
    if(sum(ind)==0)
        break
    end
    starti(ind)=starti(ind)+1;
end

if(cutoffduration>1)
    starti=starti(duration>cutoffduration);
    endi=endi(duration>cutoffduration);
    duration=duration(duration>cutoffduration);
end

%Remove storms where the window of interest extends past existing data
while(starti(1)-maxwidth<1)
    starti(1)=[]; endi(1)=[];
end
while(starti(end)+maxwidth>=length(FILLED(:,1)))
    starti(end)=[]; endi(end)=[];
end

%Cut down found storms to only include pre-noon conditions
if(cutconditions) 
   endi(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
   duration(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
   starti(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
end