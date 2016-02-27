
%Proof of Median of medians vs block median. Takahashi 2010 Figure 11
if(MakePaperPlots && stormcase==10)
    h=figure('Visible',visible);
    plot(xa(2:end),AVMDblock)
    hold on; plot(xa(2:end),AVMDs(2:end),'r')
    plot(xa(2:end),AVMD2(2:end),'k')
    xlabel('Time from minimum of event (day)')
    ylabel('\rho_{eq} (amu/cm^3)')
    legend('Block median','MoM storm, day','MoM day, storm')
    print -depsc2 -r200 paperfigures/blockmedian.eps
    if(strcmp(visible,'off')),close(h);end;
end

%Tak2006 Fig 10 Kp vs M
if(MakePaperPlots && stormcase==1)
   h=figure('Visible',visible);
   subplot(311)
   x=starti;y=starti;
   for i=1:length(starti)
    x(i)=nanmedian(FILLED(starti(i)-3:starti(i),13) );
    y(i)=nanmean(FILLED(starti(i)-3:starti(i),32).*exp(-(3:-1:0)./3)');
   end
   plot(x,y,'.'); %KP vs M
   bins=unique(x); ymid=bins; ysd=bins;
   for i=1:length(bins)
       ymid(i)=nanmedian(y(x==bins(i)));
       ysd(i)=nanstd(y(x==bins(i)));
   end
   hold on; plot(bins, ymid,'r-.','MarkerSize',10); plot(bins, ymid-ysd,'k-.','MarkerSize',10); plot(bins, ymid+ysd,'k-.','MarkerSize',10);
   title(sprintf('%d events from GOES-%d, %d-%d',length(~isnan(y)),satnum,sy,ey))
   legend('3 hour means')
   
   subplot(312) %Go from 1 hour to 1 day
   for i=1:length(starti)
    x(i)=nanmedian(FILLED(starti(i)-24:starti(i),13));
    y(i)=nanmean(FILLED(starti(i)-24:starti(i),32).*exp(-(24:-1:0)./24)');
   end
   plot(x,y,'.');
   bins=unique(x); ymid=bins; ysd=bins;
   for i=1:length(bins)
       ymid(i)=nanmedian(y(x==bins(i)));
       ysd(i)=nanstd(y(x==bins(i)));
   end
   hold on; plot(bins, ymid,'r-.','MarkerSize',10); plot(bins, ymid-ysd,'k-.','MarkerSize',10); plot(bins, ymid+ysd,'k-.','MarkerSize',10);
   legend('1 day means')
   
   subplot(313) %1 hour to 3 days
   for i=1:length(starti)
    x(i)=nanmedian(FILLED(starti(i)-24*3:starti(i),13));
    y(i)=nanmean(FILLED(starti(i)-24*3:starti(i),32).*exp(-(24*3:-1:0)./(24*3))');
   end
   plot(x,y,'.');
   bins=unique(x); ymid=bins; ysd=bins;
   for i=1:length(bins)
       ymid(i)=nanmedian(y(x==bins(i)));
       ysd(i)=nanstd(y(x==bins(i)));
   end
   hold on; plot(bins, ymid,'r-.','MarkerSize',10); plot(bins, ymid-ysd,'k-.','MarkerSize',10); plot(bins, ymid+ysd,'k-.','MarkerSize',10);
   xlabel('Kp')
   ylabel('M (amu)')
   legend('3 day means')
   
   print -depsc2 -r200 paperfigures/KpvsM.eps
   if(strcmp(visible,'off')),close(h);end;


%%%%%%%%%%%
%Tak 2010 Fig 13

   h=figure('Visible',visible);
    NewTime=FILLEDTime(1):24*27*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    x=interptest(FILLEDTime,FILLED(:,30),NewTime);
    y=log10(interptest(FILLEDTime,MassDensitySpline',NewTime));
    cc=corrcoef(x,y,'rows','pairwise');
    [AX,H1,H2]=plotyy(NewTime,x,NewTime,y,'plot','plot');
    set(H1,'marker','.','color','blue'); set(AX(1),'YColor','r'); set(AX(2),'XTick',[]);
    set(H2,'marker','.','color','red'); set(AX(2),'YColor','b');
    ylim(AX(1),[0,300])
    ylim(AX(2),[0.5,1.5])
    ylabel(AX(1),'F_{10.7} (s.f.u.)','FontSize',BigFont); ylabel(AX(2),'GOES-6 log_{10}[\rho_{eq} (amu/cm^2)]','FontSize',BigFont);
    set(AX(1),'YTick',0:50:300);
    %set(AX(2),'YTick',.5:0.25:1.5);
    xlabel('Year','FontSize',BigFont);
    title('27-day medians','FontSize',BigFont)
    linkaxes(AX,'x')
    datetick('keeplimits');
    grid on
    print('-depsc2', '-r200', sprintf('paperfigures/F107MD27d-GOES%d.eps',satnum))
    print('-dpdf', sprintf('paperfigures/F107MD27d-GOES%d.pdf',satnum))
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/F107MD27d-GOES%d.png',satnum))
    if(strcmp(visible,'off')),close(h);end;
    
    h=figure('Visible',visible);
    NewTime=FILLEDTime(1):24*(FILLEDTime(2)-FILLEDTime(1)):FILLEDTime(end);
    x=interptest(FILLEDTime,FILLED(:,30),NewTime);
    y=log10(interptest(FILLEDTime,MassDensitySpline',NewTime));
    cc=corrcoef(x,y,'rows','pairwise');
    [AX,H1,H2]=plotyy(NewTime,x,NewTime,y,'plot','plot');
    set(H1,'marker','.','color','red'); set(AX(1),'YColor','r'); set(AX(2),'XTick',[]);
    ylim(AX(2),[0.5,1.5])
    ylabel(AX(1),'F_{10.7\_1d}','FontSize',BigFont); ylabel(AX(2),'log(\rho_{eq\_1d})','FontSize',BigFont);
    xlabel('Year','FontSize',BigFont);
    title(sprintf('Data comparison for GOES %d - CC: %2.2f',satnum,cc(1,2)),'FontSize',BigFont)
    linkaxes(AX,'x')
    datetick('keeplimits');
    grid on
    print('-depsc2', '-r200', sprintf('paperfigures/F107MD1d-GOES%d.eps',satnum))
    print('-dpng', '-r200', sprintf('paperfigures/PNGs/F107MD1d-GOES%d.png',satnum))
    if(strcmp(visible,'off')),close(h);end;
end