function [s, st]=twobinplot(plotter,sorter,events,timewidth,LongTimeScale,plotthresh,names,units,timerange,visible)

if nargin<4
    fprintf('Usage:')
    fprintf('h=twobinplot(plotter,sorter,events,timewidth,LongTimeScale,plotthresh,names,units,timerange,visible)')
    fprintf('names={plotvar name; sortvar name; event name}')
    return
end

plotname=char(names(1));
sortname=char(names(2));
eventname=char(names(3));
plotunits=char(units(1));
sortunits=char(units(2));
eventunits=char(units(3));
safeplotname=regexprep(plotname,'[^a-zA-Z0-9]','');
safesortname=regexprep(sortname,'[^a-zA-Z0-9]','');
safeeventname=regexprep(eventname,'[^a-zA-Z0-9]','');
safethresh=regexprep(plotthresh,'[^a-zA-Z0-9]','');

xa=(-timewidth:LongTimeScale:timewidth*2);

medsplit=nanmedian(sorter(events));

HighIndex=sorter(events)>medsplit;
LowIndex=sorter(events)<=medsplit;
HighIndices=HighIndex'.*(1:length(HighIndex));
LowIndices=LowIndex'.*(1:length(LowIndex));
HighIndices(HighIndices==0)=[];
LowIndices(LowIndices==0)=[];

high=nanmedian(plotter(HighIndex,:),1);
low=nanmedian(plotter(LowIndex,:),1);

highbase=nanmedian([high(xa<-6) high(xa>12)]);
lowbase=nanmedian([low(xa<-6) low(xa>12)]);

for i=1:100
    bootshigh(i,:)=nanmedian(plotter(randsample(HighIndices,size(plotter,1),true),:));
    bootslow(i,:)=nanmedian(plotter(randsample(LowIndices,size(plotter,1),true),:));
    %bootstdhigh(i,:)=nanstd(plotter(randsample(HighIndex,size(plotter,1),true),:));
    %bootstdlow(i,:)=nanstd(plotter(randsample(LowIndex,size(plotter,1),true),:));
    ss(i)=sum((bootshigh(i,:)-highbase)-(bootslow(i,:)-lowbase));
    sts(i)=std((bootshigh(i,:)-highbase)-(bootslow(i,:)-lowbase));
end

s=nanmean(ss);
st=nanmean(sts);


%s=sum((high-highbase)-(low-lowbase));
%st=std((high-highbase)-(low-lowbase));

HighMatBars=[nanmedian(plotter(HighIndex,:))-highbase-nanstd(plotter(HighIndex,:))./sqrt(sum(~isnan(plotter(HighIndex,:)))) ; nanmedian(plotter(HighIndex,:))-highbase+nanstd(plotter(HighIndex,:))./sqrt(sum(~isnan(plotter(HighIndex,:))))];
LowMatBars=[nanmedian(plotter(LowIndex,:))-lowbase-nanstd(plotter(LowIndex,:))./sqrt(sum(~isnan(plotter(LowIndex,:)))) ; nanmedian(plotter(LowIndex,:))-lowbase+nanstd(plotter(LowIndex,:))./sqrt(sum(~isnan(plotter(LowIndex,:))))];

%Bootstrap instead

%bootstdhigh=nanmedian(bootstdhigh);
%bootstdlow=nanmedian(bootstdlow);
%HighMatBars=[nanmedian(plotter(HighIndex,:))-highbase-bootstdhigh ; nanmedian(plotter(HighIndex,:))-highbase+bootstdhigh];
%LowMatBars=[nanmedian(plotter(LowIndex,:))-lowbase-bootstdlow ; nanmedian(plotter(LowIndex,:))-lowbase+bootstdlow];

h=figure('visible',visible);
plot(xa,high-highbase,'b'); hold on; plot(xa,low-lowbase,'r');
plot(xa,HighMatBars,'b-.'); plot(xa,LowMatBars,'r-.');
ylabel(sprintf('%s (%s)',plotname,plotunits))
xlabel('Time from start of event (hour)')
set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
legend(sprintf('%s>%2.0f',sortname,medsplit),sprintf('%2.0f>%s',medsplit,sortname),'Location','NorthEast');
title(sprintf('%d evenly-binned events of %s %s (%s) for %d-%d, baseline removed',length(events),eventname,plotthresh,eventunits,timerange(1),timerange(2)))
print('-depsc2','-r200',sprintf('figures/HighLowBase%s%s-%s%s.eps',safesortname,safeplotname,safeeventname,safethresh))
print('-dpng','-r200',sprintf('figures/HighLowBase%s%s-%s%s.png',safesortname,safeplotname,safeeventname,safethresh))

%{
if(sum(sum(isnan(plotter)))) %If any nan points
    h=figure('Visible',visible);
    plot(xa,sum(~isnan(plotter(sorter(events)>highsplit,:))))
    hold on;
    plot(xa,sum(~isnan(plotter(sorter(events)>medsplit & sorter(events)<highsplit,:))),'g')
    plot(xa,sum(~isnan(plotter(sorter(events)<medsplit & sorter(events)>lowsplit,:))),'r')
    plot(xa,sum(~isnan(plotter(sorter(events)<lowsplit,:))),'c')
    legend(sprintf('%s>%2.0f',sortname,highsplit),sprintf('%2.0f>%s>%2.0f',highsplit,sortname,medsplit),sprintf('%2.0f>%s>%2.0f',medsplit,sortname,lowsplit),sprintf('%2.0f>%s',lowsplit,sortname),'Location','NorthEast');
    xlabel('Time from start of event (hour)')
    ylabel('Valid hourly data points')
    title(sprintf('%d evenly-binned events of %s %s (%s) for %d-%d',length(events),eventname,plotthresh,eventunits,timerange(1),timerange(2)))
    axis tight;
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
    print('-depsc2','-r200',sprintf('paperfigures/HighLow%s%s-%s%s-valid.eps',safesortname,safeplotname,safeeventname,safethresh))
    print('-dpng','-r200',sprintf('paperfigures/PNGs/HighLow%s%s-%s%s-valid.png',safesortname,safeplotname,safeeventname,safethresh))
end
%}
