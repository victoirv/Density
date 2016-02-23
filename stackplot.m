function stackplot(Time,Data, labels, satnum, ylims, datelimits)
%To ignore ylimits, pass ylims=0
%TODO: include filename argument?

if(nargin < 5)
    ylims=0;
    datelimits=[0 datenum('Jan-01-3001')];
elseif nargin < 6
    datelimits=[0 datenum('Jan-01-3001')];
end

h=figure;%('Visible',visible);
orient tall;
hold on;
numvars=length(labels);

Data(Time<datelimits(1),:)=[];
Time(Time<datelimits(1))=[]; 
Data(Time>datelimits(2),:)=[];
Time(Time>datelimits(2))=[];

for i=1:numvars
    h(i)=subplot('position',subplotstack(numvars,i));
    plot(Time,Data(:,i),'.','MarkerSize',12); 
    text(0.01,0.85,labels{i},'Units','normalized','FontSize',14); 
    if(ylims~=0)
        ylim(ylims(i,:));
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
    xlabel('Time')
    
    sy=year(Time(1));
    ey=year(Time(end));
    print('-depsc2','-r200',sprintf('paperfigures/alldata-GOES%d-%d-%d.eps',satnum,sy,ey));
    print('-dpng','-r200',sprintf('paperfigures/PNGs/alldata-GOES%d-%d-%d.png',satnum,sy,ey));