function [starti,endi,duration]=FindStorms(storms,cutoffduration,cutconditions,maxwidth,MLTFit,FILLED,DstCut)
%If reproduction of DST shifted events is needed, pass FILLED

starti=find(storms>0);
endi=find(storms<0)-1;
duration=endi-starti+1;

%Shift event start points to next local DST minimum 
while(1 && cutconditions && DstCut~=0) 
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
while(starti(end)+maxwidth>=length(storms)-1)
    starti(end)=[]; endi(end)=[];
end

tooclose=find(diff(starti)<23)+1;
starti(tooclose)=[]; %remove storms within 10 hours of each other
endi(tooclose)=[];
duration(tooclose)=[];

%Cut down found storms to only include pre-noon conditions
if(cutconditions) 
    
    
   endi(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
   duration(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
   starti(MLTFit(starti)>12 | MLTFit(starti)<6)=[];
   
   drop=[];
   for i=1:length(starti)
       if(sum(FILLED([starti(i)-20 : starti(i)-1],15) < FILLED(starti(i),15))>0)
           drop=[drop i];
       end
   end
   starti(drop)=[];
   endi(drop)=[];
   duration(drop)=[];
   
end