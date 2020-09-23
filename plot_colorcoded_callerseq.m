clear all;
close all;
colors=[146 192 215;183 251 254;6 127 122;44 100 87;0 40 22]./255;
colors2=[188 134 150;254 200 126;254 235 167;126 214 190;85 124 155]./255;
colors3=[92 48 107;195 49 122;211 200 170;200 176 6;73 165 154]./255;
colors4=[197 5 93;249 118 46;248 182 36;255 237 215;172 232 242]./255;
colors5=[167 124 167; 70 37 90;87 122 164;245 235 210;192 168 144]./255;
load('C:\Users\Daria\audio_stuff\all_cut_call_data.mat')
icount=0;
for i=1:length(all_cut_call_data)
    if length(all_cut_call_data(i).batNum)==1
        if str2num(cell2mat(all_cut_call_data(i).batNum))>1000
            icount=icount+1;
            batids(icount)=str2num(cell2mat(all_cut_call_data(i).batNum));
            alltreat{icount}=all_cut_call_data(i).treatment;
            dat=all_cut_call_data(i).fName;
            if size(dat,2)<50
                dat=dat{1};
            end
            dat=str2num(dat(42:49));
            all(icount)=dat;
            pos(icount,:)=(all_cut_call_data(i).corrected_callpos);
            
            
        end
    end
    
end
expdays=(unique(all));
exd=zeros(length(all),1);
for i=1:length(expdays)
    indx=find(all==expdays(i));
    exd(indx)=i;
end
batcolors=ones(length(batids),3).*.7;
trts=zeros(length(batids),1);
for run=1:2
    if run==2
        %deaf
        IDs=[14463;71043;71047;71351;71354];
        val=colors2;%[0 0 1;1 0 0;0 1 0;1 1 0;1 0 1];
    else
        %saline
        IDs=[11648;14461;14464;65696;71353];
        val=colors5;%[0 0 1;1 0 0;0 1 0;1 1 0;1 0 1];
    end
    for bn=1:length(IDs)
        indx=find(batids==IDs(bn));
        batcolors(indx,1)=val(bn,1);
        batcolors(indx,2)=val(bn,2);
        batcolors(indx,3)=val(bn,3);
        trts(indx)=run;
    end
end
figure;
for n=1:size(batcolors,1)
    trans=0.1;
    
    if n>1
        dpos(n-1)=abs(diff([pos(n,1) pos(n-1,2)]));
        if abs(diff([pos(n,1) pos(n-1,2)]))<500
            trans=1;
        elseif abs(diff([pos(n-1,2) pos(n,1)]))>=500 && abs(diff([pos(n-1,2) pos(n,1)]))<3e3
            trans=.5;
        end
    end
    [fillhandle,msg]=jbfill((n-1):n,[0 0]+(trts(n)-1),[1 1]+(trts(n)-1),batcolors((n),:),batcolors((n),:),trans,trans);
    hold on;
end
for n=1:length(expdays)
    indx=find(exd==n);
    plot([indx(end) indx(end)],[-1 2],'k:','linewidth',2);
    hold on;
end
set(gca,'ylim',[-1 2],'ytick',(-.5:1:1.5),'yticklabel',{'other';'saline';'deaf'},'xlim',[0 indx(end)])
xlabel('Call # (dashed lines = session separation)')

bats=unique(batids);
for bn=1:length(bats)
    indx=find(batids==bats(bn));
    for neighb=1:length(bats)
        n=0;
        for i=1:length(indx)
            if (indx(i)+1)<=length(batids)
                if batids(indx(i)+1)==bats(neighb)
                    n=n+1;
                end
            end
        end
        neighbc(bn,neighb,1)=n;
        n=0;
        for i=1:length(indx)
            if (indx(i)-1)>=1
                if batids(indx(i)-1)==bats(neighb)
                    n=n+1;
                end
            end
        end
        neighbc(bn,neighb,2)=n;
    end
end
%deaf
IDs=[14463;71043;71047;71351;71354];
for pl1=1:size(neighbc,1)
    for pl2=1:size(neighbc,2)
        if squeeze(neighbc(pl1,pl2,1))~=0
            if ismember(bats(pl1),IDs)
              figure(10);  plot(pl1,pl2,'r.','Markersize',squeeze(neighbc(pl1,pl2,1))/20)
            else
              figure(10);  plot(pl1,pl2,'k.','Markersize',squeeze(neighbc(pl1,pl2,1))/20)
+            end
            hold on;
        end
        if squeeze(neighbc(pl1,pl2,2))~=0
            if ismember(bats(pl2),IDs)
             figure(11);   plot(pl1,pl2,'r.','Markersize',squeeze(neighbc(pl1,pl2,2))/20)
            else
             figure(11);   plot(pl1,pl2,'k.','Markersize',squeeze(neighbc(pl1,pl2,2))/20)
            end
            hold on;
        end
        
    end
end
figure(10);
set(gca,'xlim',[0 length(bats)+1],'xtick',(1:length(bats)),'xticklabel',bats,...
    'ylim',[0 length(bats)+1],'ytick',(1:length(bats)),'yticklabel',bats)
xtickangle(90)
xlabel('caller')
ylabel('following caller')
figure(11);
set(gca,'xlim',[0 length(bats)+1],'xtick',(1:length(bats)),'xticklabel',bats,...
    'ylim',[0 length(bats)+1],'ytick',(1:length(bats)),'yticklabel',bats)
xtickangle(90)
xlabel('following caller')
ylabel('caller')

