clear
addpath('~/git/m-rsw/time');
addpath('~/git/m-rsw/stats');

load 1dData.mat
Nv = NUsed;
Nd = 4;
Nr = floor(length(F107d)/Nd);

F107 = F107d(1:Nd*Nr);
Dst  = Dstd(1:Nd*Nr);
MD    = MDd(1:Nd*Nr);
Nv    = Nv(1:Nd*Nr);

mode = 'none';

if strcmp(mode,'none')
    MDdr   = nanmedian(reshape(MD,Nd,Nr),1)';
    F107r  = nanmedian(reshape(F107,Nd,Nr),1)';
    Dstr   = nanmedian(reshape(Dst,Nd,Nr),1)';
    Dstr   = nanstd(reshape(F107,Nd,Nr),1)';
    ls = {'\rho_{eq}','F10.7','D_{st}','N_{valid}'};
end

I    = find(isnan(MDdr) == 0);
MDdi = interp1(I,MDdr(I),1:length(MDdr))';
MDdi(isnan(MDdi)) = 0;

corrcoef(log10(MDdr),F107r,'rows','pairwise')

[T,X] = time_delay(log10(MDdr),F107r,10,-1,0,'pad');
LIN   = basic_linear(X,T)
MDp   = basic_linear(X,LIN.Weights,'predict');
MDe   = log10(MDdr(2:end))-MDp;
MDe(isnan(MDe)) = 0;
figure(1);clf;grid on;hold on;
plot(MDe,'LineWidth',2);
plot(-(Dstr)/100,'LineWidth',1);

legend('Error','Dst')
[T,X] = time_delay(MDe,Dstr(2:end),10,-1,0,'pad');
LIN   = basic_linear(X,T)
corrcoef(abs(MDe),Dstr(2:end))

break
[c,lags] = xcorr(MDe,'unbiased');
plot(lags,c,'.');
break

LIN  = basic_linear(F107t,log10(MDt))
LIN  = basic_linear(F107m,log10(MDm))
MDdp = basic_linear(F107m,LIN.Weights,'predict');
MDe  = MDm-MDdp;

if (1)
figure(5);clf;hold on;grid on;box on;
   t = [0:length(F107)-1]';
   plot(t,MDe)
   plot(t,100+F107-F107m,'k','LineWidth',2);  
end

MDd = MDe;
I    = find(isnan(MDd) == 0);
MDdi = interp1(I,MDd(I),1:length(MDd))';
MDdi(isnan(MDdi)) = 0;

[I_F107d,f] = periodogramraw(F107d,'fast');
I_Dstd      = periodogramraw(Dstd,'fast');
I_MDdi      = periodogramraw(MDdi,'fast');
I_Nv        = periodogramraw(Nv,'fast');
lw          = {'LineWidth',2};

W = [1:2*Nd];
for w = 1:length(W);
    S_F107d = block_std(F107d,w);
    S_Dstd  = block_std(Dstd,w);
    S_MDdi  = block_std(MDdi,w);
    S_Nv    = block_std(Nv,w);
    
    M_F107d = block_mean(F107d,w);
    M_Dstd  = block_mean(Dstd,w);
    M_MDdi  = block_mean(MDdi,w);
    M_Nv    = block_mean(Nv,w);    

    tmp = corrcoef(S_Nv,S_MDdi);
    S_ccNv(w) = tmp(2);
    tmp = corrcoef(S_F107d,S_MDdi);
    S_ccF107(w) = tmp(2);
    tmp = corrcoef(S_Dstd,S_MDdi);
    S_ccDst(w) = tmp(2);

    tmp = corrcoef(M_Nv,M_MDdi);
    M_ccNv(w) = tmp(2);
    tmp = corrcoef(M_F107d,M_MDdi);
    M_ccF107(w) = tmp(2);
    tmp = corrcoef(M_Dstd,M_MDdi);
    M_ccDst(w) = tmp(2);

    %fprintf('%f %.2f %.2f\n',w,cc1(w),cc2(w));
end

plot(S_F107d,S_MDdi,'.','MarkerSize',20);

if (0)
Nl = 10
[T,X] = time_delay(MDd,F107d,Nl,-1);
LIN = basic_linear(X,T)
[T,X] = time_delay(MDd,Dstd,Nl,-1);
LIN = basic_linear(X,T)
[T,X] = time_delay(MDd,Nv,Nl,-1);
LIN = basic_linear(X,T)
end

Nl = 1
[T,X] = time_delay(S_MDdi,[S_F107d],Nl,-1);
LIN = basic_linear(X,T)
[T,X] = time_delay(S_MDdi,[S_F107d,S_Dstd],Nl,-1);
LIN = basic_linear(X,T)

t = [0:length(F107d)-1]/Nd;  
figure(1);clf;grid on;hold on;box on;
    plot(t,0+MDdi/max(F107d),'b',lw{:});
    plot(t,1+F107d/max(F107d),'k',lw{:});
    plot(t,2+Dstd/max(F107d),'r',lw{:});
    plot(t,3+Nv/max(Nv),'g',lw{:});
    xlabel(sprintf('time [%d day increments]',Nd));
    set(gca,'XTick',[0:5:115])
    if (0)
    xl = get(gca,'XTickLabel');
    for i = 1:length(xl)
	if (mod(i,5) ~= 0)
	    xl{i} = '';
	end
    end
    set(gca,'XTickLabel',xl);
    end
    for i = 1:length(ls)
	text(t(end),i-1,sprintf(' %s',ls{i}));
    end
    set(gca,'XLim',[0,t(end)]);

figure(2);clf;hold on;grid on;box on;
    plot(1./f(2:end),0+I_MDdi(2:end)/max(I_MDdi(2:end)),'b',lw{:});
    plot(1./f(2:end),1+I_F107d(2:end)/max(I_F107d(2:end)),'k',lw{:});
    plot(1./f(2:end),2+I_Dstd(2:end)/max(I_Dstd(2:end)),'r',lw{:});
    plot(1./f(2:end),3+I_Nv(2:end)/(max(I_Nv(2:end))),'g',lw{:});
    xlabel('Period [days]');
    %set(gca,'XLim',[0,100]);
    %set(gca,'XTick',[2,[3:3:30]])
    %set(gca,'XLim',[2,30]);
    set(gca,'YTick',[]);
    %xlabel('f [1/days]');
    title('Raw Periodogram (arb. units)');
    for i = 1:length(ls)
	text(30,i-1,sprintf(' %s',ls{i}));
    end

figure(3);;clf;hold on;grid on;box off;
  plot(W,M_ccNv,'g','LineWidth',2);
  text(W(end),M_ccNv(end),sprintf('cc(%s,%s)',ls{4},ls{1}))
  plot(W,M_ccDst,'r','LineWidth',2); 
  text(W(end),M_ccDst(end),sprintf('cc(%s,%s)',ls{3},ls{1}))
  plot(W,M_ccF107,'k','LineWidth',2);
  text(W(end),M_ccF107(end),sprintf('cc(%s,%s)',ls{2},ls{1}))
  xlabel('Window width [days]');
  %ylabel('cc','Rotation',0,'HorizontalAlignment','Right');
    
figure(4);;clf;hold on;grid on;box off;
  plot(W,S_ccNv,'g','LineWidth',2);
  text(W(end),S_ccNv(end),sprintf('cc(\\sigma_{%s},%s)',ls{4},ls{1}))
  plot(W,S_ccDst,'r','LineWidth',2); 
  text(W(end),S_ccDst(end),sprintf('cc(\\sigma_{%s},%s)',ls{3},ls{1}))
  plot(W,S_ccF107,'k','LineWidth',2);
  text(W(end),S_ccF107(end),sprintf('cc(\\sigma_{%s},%s)',ls{2},ls{1}))
  xlabel('Window width [days]');
  %ylabel('cc','Rotation',0,'HorizontalAlignment','Right');

break  
for i = 1:37*3
    [tmp,lags] = xcorr(F107dr(:,i),MDdir(:,i),27*1);
    C_F107d_MDdi(:,i) = tmp;
end

[C,lags] = xcorr(F107dr(:,i),MDdir(:,i),27*1);
figure(3);clf;grid on;hold on;
    plot(lags,mean(C_F107d_MDdi,2))

figure(4)
    plot(LIN.Weights)