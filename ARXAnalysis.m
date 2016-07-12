function ARXAnalysis(FILLED,MD,headers,satnum)

loops=100;
trainpercent=0.50;

table=fopen(sprintf('tables/ARXtable-GOES%d.txt',satnum),'w');
fprintf(table,'<pre>\n');
fprintf(table,'Vars \t \t  CC1  CC24 CC24+1\n');

MD=MD(:); %Make column vector

FILLED=FILLED(:,[1:15 29 30 31]);
headers=headers([1:15 29 30 31]);

for permlen=1:2
    Perms=combnk(1:length(headers),permlen);
    for i=1:length(Perms)
        x=FILLED(:,Perms(i,:));
        x=reshape(x,length(x),[]);
        y=MD;
        Z=[x y];
        for col=1:min(size(Z))
            Z(isnan(Z(:,col)),:)=[];
        end
        x=Z(:,1:end-1);
        y=Z(:,end);

        [~,~,~,~, corr1, eff]=IR(y,x,0,1,0,0,loops);
        [~,~,~,~, corr12, eff]=IR(y,x,0,24,0,0,loops);
        [~,~,~,~, corr241, eff]=IR(y,x,1,24,0,0,loops);
        
        
        
        combohead=sprintf('%s',headers{Perms(i,1)});
        if(length(Perms(i,:))>1)
            for hc=2:length(Perms(i,:))
                combohead=sprintf('%s+%s',combohead,headers{Perms(i,hc)});
            end
        end
        fprintf(table,'%s      \t- %+2.2f %+2.2f %+2.2f\n',combohead,corr1,corr12,corr241);
    end
end


Perms=combnk(1:length(headers),length(headers));
        x=FILLED(:,Perms(1,:));
        x=reshape(x,length(x),[]);
        y=MD;
        Z=[x y];
        for col=1:min(size(Z))
            Z(isnan(Z(:,col)),:)=[];
        end
        x=Z(:,1:end-1);
        y=Z(:,end);

        [~,~,~,~, corr1, eff]=IR(y,x,0,1,0,0,loops);
        [~,~,~,~, corr12, eff]=IR(y,x,0,24,0,0,loops);
        [~,~,~,~, corr241, eff]=IR(y,x,1,24,0,0,loops);
fprintf(table,'All\t- %+2.2f %+2.2f %+2.2f\n',corr1,corr12,corr241);

fclose(table);