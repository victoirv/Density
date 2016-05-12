function h=binplot(plotter,sorter,events,timewidth,LongTimeScale,plotthresh,names,units,timerange,satnum,visible)

if nargin<4
    fprintf('Usage:')
    fprintf('h=binplot(plotter,sorter,events,timewidth,LongTimeScale,names)')
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
highsplit=nanmedian(sorter(events(sorter(events)>medsplit)));
lowsplit=nanmedian(sorter(events(sorter(events)<medsplit)));

midhigh=nanmedian(plotter(sorter(events)>medsplit & sorter(events)<highsplit,:),1);
midlow=nanmedian(plotter(sorter(events)<medsplit & sorter(events)>lowsplit,:),1);
high=nanmedian(plotter(sorter(events)>highsplit,:),1);
low=nanmedian(plotter(sorter(events)<lowsplit,:),1);
HighIndex=sorter(events)>highsplit;
LowIndex=sorter(events)<lowsplit;
HighMatBars=[nanmedian(plotter(HighIndex,:))-nanstd(plotter(HighIndex,:))./sqrt(sum(~isnan(plotter(HighIndex,:)))) ; nanmedian(plotter(HighIndex,:))+nanstd(plotter(HighIndex,:))./sqrt(sum(~isnan(plotter(HighIndex,:))))];
LowMatBars=[nanmedian(plotter(LowIndex,:))-nanstd(plotter(LowIndex,:))./sqrt(sum(~isnan(plotter(LowIndex,:)))) ; nanmedian(plotter(LowIndex,:))+nanstd(plotter(LowIndex,:))./sqrt(sum(~isnan(plotter(LowIndex,:))))];

h=figure('visible',visible);
plot(xa,[high; midhigh; midlow; low]);
hold on; plot(xa,HighMatBars,'b-.'); plot(xa,LowMatBars,'c-.');
ylabel(sprintf('%s (%s)',plotname,plotunits))
xlabel('Time from start of event (hour)')
set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
[~,objh]=legend(sprintf('%s>%2.0f',sortname,highsplit),sprintf('%2.0f>%s>%2.0f',highsplit,sortname,medsplit),sprintf('%2.0f>%s>%2.0f',medsplit,sortname,lowsplit),sprintf('%2.0f>%s',lowsplit,sortname),'Location','NorthEast');
set(objh,'linewidth',2);
title(sprintf('%d evenly-binned events of %s %s (%s) for GOES%d: %d-%d',length(events),eventname,plotthresh,eventunits,satnum,timerange(1),timerange(2)))
print('-depsc2','-r200',sprintf('figures/HighLow%s%s-%s%s-GOES%d-%d-%d.eps',safesortname,safeplotname,safeeventname,safethresh,satnum,timerange(1),timerange(2)))
print('-dpng','-r200',sprintf('figures/PNGs/HighLow%s%s-%s%s-GOES%d-%d-%d.png',safesortname,safeplotname,safeeventname,safethresh,satnum,timerange(1),timerange(2)))
if(strcmp(visible,'off')),close(h);end;

if(sum(sum(isnan(plotter)))) %If any nan points
    h=figure('Visible',visible);
    plot(xa,sum(~isnan(plotter(sorter(events)>highsplit,:))))
    hold on;
    plot(xa,sum(~isnan(plotter(sorter(events)>medsplit & sorter(events)<highsplit,:))),'g')
    plot(xa,sum(~isnan(plotter(sorter(events)<medsplit & sorter(events)>lowsplit,:))),'r')
    plot(xa,sum(~isnan(plotter(sorter(events)<lowsplit,:))),'c')
    [~,objh]=legend(sprintf('%s>%2.0f',sortname,highsplit),sprintf('%2.0f>%s>%2.0f',highsplit,sortname,medsplit),sprintf('%2.0f>%s>%2.0f',medsplit,sortname,lowsplit),sprintf('%2.0f>%s',lowsplit,sortname),'Location','NorthEast');
    set(objh,'linewidth',2);
    xlabel('Time from start of event (hour)')
    ylabel('Valid hourly data points')
    title(sprintf('%d evenly-binned events of %s %s (%s) for GOES%d: %d-%d',length(events),eventname,plotthresh,eventunits,satnum,timerange(1),timerange(2)))
    axis tight;
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
    print('-depsc2','-r200',sprintf('figures/HighLow%s%s-%s%s-GOES%d-%d-%d-valid.eps',safesortname,safeplotname,safeeventname,safethresh,satnum,timerange(1),timerange(2)))
    print('-dpng','-r200',sprintf('figures/PNGs/HighLow%s%s-%s%s-GOES%d-%d-%d-valid.png',safesortname,safeplotname,safeeventname,safethresh,satnum,timerange(1),timerange(2)))
    if(strcmp(visible,'off')),close(h);end;
end
