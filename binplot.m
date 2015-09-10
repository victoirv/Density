function h=binplot(plotter,sorter,events,timewidth,LongTimeScale,plotthresh,names,units,visible)

if nargin<4
    fprintf('Usage:')
    fprintf('h=binplot(plotter,sorter,events,timewidth,LongTimeScale,names)')
    fprintf('names=[plotvar name; sortvar name]')
    return
end
    plotname=char(names(1));
    sortname=char(names(2));
    plotunits=char(units(1));
    sortunits=char(units(2));
    safeplotname=regexprep(plotname,'[^a-zA-Z0-9]','');
    safesortname=regexprep(sortname,'[^a-zA-Z0-9]','');
    
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
    ylabel(sprintf('%s %s',plotname,plotunits))
    xlabel('Time from start of event (hour)')
    set(findobj('type','axes'),'xgrid','on','ygrid','on','box','on','xtick',[-timewidth:timewidth/2:timewidth*2]./LongTimeScale)
    legend(sprintf('%s>%2.0f',sortname,highsplit),sprintf('%2.0f>%s>%2.0f',sortname,highsplit,medsplit),sprintf('%2.0f>%s>%2.0f',sortname,medsplit,lowsplit),sprintf('%2.0f>%s',sortname,lowsplit),'Location','NorthEast');
    title(sprintf('%d evenly-binned events of %s < %d %s for %d-%d',length(events),plotname,plotthresh,plotunits,year(OMNITime(1)),year(OMNITime(end))))
    print('-depsc2','-r200','paperfigures/HighLow%s%s.eps',safesortname,safeplotname)
    print('-dpng','-r200','paperfigures/PNGs/HighLow%s%s.png',safesortname,safeplotname)