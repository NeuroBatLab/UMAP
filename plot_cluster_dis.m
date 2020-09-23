function plot_hopkins_dis(pathname,bin,compN)

load([pathname(1:end-4) '_EVAL.mat'])
cluscomp=ClusComps{compN};
for rp=1:size(cluscomp,2)
    fh=figure('Name',['Cluster dist.-Comp#' num2str(compN) '-Rep#' num2str(rp)],'NumberTitle','off');
    count=0;
    for bnr=1:size(cluscomp,1)
        alli=[];
        for mi=1:2
            indx=find(~isnan(cluscomp(bnr,rp,mi,:)));
            if ~isempty(indx)
                alli(mi)=indx(end);
            end
        end
        if ~isempty(alli)
            count=count+1;
            figure(fh);subplot(ceil(size(cluscomp,1)/ceil(size(cluscomp,1)/2)),ceil(size(cluscomp,1)/2),count);
            bar(1:max(alli),squeeze(cluscomp(bnr,rp,:,1:max(alli)))')
            set(gca,'ylim',[0 1])
            title(['Bin: ' num2str(bin(bnr))])
            if count==1
                legend('Control','Deaf')
            end
            if count==4
                ylabel('Cluster (%)')
            elseif count==5
                xlabel('Cluster #')
            end
        end
    end
end