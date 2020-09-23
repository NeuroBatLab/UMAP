function plot_init_umap_result(pathname,PARAS,bin)

fh2=figure('Name','Hopkins results','NumberTitle','off');
fh3=figure('Name','dbScan results','NumberTitle','off');
fh4=figure('Name','hdbScan results','NumberTitle','off');
load([pathname(1:end-4) '_EVALinit.mat'])
for bnr=1:size(hops,1)
    count=0;
    for md=1:size(hops,3)
        for nn=1:size(hops,4)
            count=count+1;
            figure(fh2);
            subplot(ceil(size(hops,1)/ceil(size(hops,1)/2)),ceil(size(hops,1)/2),bnr);hold on;
            plot(count,squeeze(hops(bnr,:,md,nn,1)),'k.','Markersize',10)
            xlab{count}=[num2str(PARAS.min_dis(md)) '\_' num2str(PARAS.n_neighbors(nn))];
            figure(fh3);
            subplot(ceil(size(clN,1)/ceil(size(clN,1)/2)),ceil(size(clN,1)/2),bnr);hold on;
            plot(count,squeeze(clN(bnr,:,md,nn,1)),'k.','Markersize',10)
            figure(fh4);
            subplot(ceil(size(clN,1)/ceil(size(clN,1)/2)),ceil(size(clN,1)/2),bnr);hold on;
            plot(count,squeeze(clN(bnr,:,md,nn,2)),'k.','Markersize',10)
        end
    end
    for pl=1:3
        if pl==1
            figure(fh2);
        elseif pl==2
            figure(fh3);
        else
            figure(fh4);
        end
        set(gca,'xlim',[0 count+1],'xtick',(1:1:count),'xticklabel',xlab)
        xtickangle(90)
        title(['Bin: ' num2str(bin(bnr))])
        if pl==1
            set(gca,'ylim',[0.4 1])
        else
            set(gca,'ylim',[0 50])
        end
        if bnr==size(hops,1)/2+1
            if pl==1
                ylabel('Hopkins value')
            elseif pl==2
                ylabel('Sum db-Cluster (#)')
            else
                ylabel('Sum hdb-Cluster (#)')
            end
        elseif bnr>size(hops,1)-ceil(size(hops,1)/2)
            xlabel('MD-NN pair')
        end
    end
end
