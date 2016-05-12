function stackplot(Time,Data, labels, satnum, ThreshLines,ylims, datelimits,visible)
%To ignore ylimits, pass ylims=0
%TODO: include filename argument?

if (nargin < 5 || isempty(ThreshLines))
    ThreshLines=[0 0]; %Since there is no zero-th plot
end
if(nargin < 6 || isempty(ylims))
    ylims=0;
end
if nargin < 7
    datelimits=[0 datenum('Jan-01-3001')]; %Arbitrary time cutoffs that will hopefully never be reached
end
if nargin < 8
    visible='off';
end



h=figure('Visible',visible);
orient tall;
hold on;
numvars=length(labels);

Data(Time<datelimits(1),:)=[];
Time(Time<datelimits(1))=[];
Data(Time>datelimits(2),:)=[];
Time(Time>datelimits(2))=[];

for i=1:numvars
    h(i)=subplot('position',subplotstack(numvars,i));
    if(i<numvars)
        plot(Time,Data(:,i),'.','MarkerSize',10);
    else
        plot(Time,Data(:,i),'r.','MarkerSize',10); %Plot last subplot in red
    end
    if(sum(ThreshLines(:,1)==i)>0)
        hold on; plot(Time,ones(1,length(Time)).*(ThreshLines(ThreshLines(:,1)==i,2)),'k-.','LineWidth',2)
    end
    text(0.01,0.85,labels{i},'Units','normalized','FontSize',14);
    if(ylims~=0)
        ylim(ylims(1,i,:));
    end
end

%hold on; plot([Time(1) Time(end)],[-40 -40],'r-.','LineWidth',4); hold off;
%hold on; plot([Time(1) Time(end)],[40 40],'b-.','LineWidth',4);

set(findobj('type','axes'),'xticklabel',{[]});
set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on')
axis tight;
datetick('x','keeplimits')
set(findobj('type','axes'),'xtick',get(h(end),'xtick'))
linkaxes(h,'x')
axis tight;
xlabel('Time','FontSize',14)

sy=year(Time(1));
ey=year(Time(end));

if(sy ~= ey)
%print('-depsc2',sprintf('paperfigures/alldata-GOES%d-%d-%d.eps',satnum,sy,ey));
print('-depsc2','-r300',sprintf('paperfigures/alldata-GOES%d-%d-%d.eps',satnum,sy,ey));
print('-dpng','-r200',sprintf('paperfigures/PNGs/alldata-GOES%d-%d-%d.png',satnum,sy,ey));
else
    sy=datestr(Time(1),'ddmmmyyyy');
    ey=datestr(Time(end),'ddmmmyyyy');
    print('-depsc2','-r300',sprintf('paperfigures/alldata-GOES%d-%s-%s.eps',satnum,sy,ey));
    print('-dpng','-r200',sprintf('paperfigures/PNGs/alldata-GOES%d-%s-%s.png',satnum,sy,ey));
end
