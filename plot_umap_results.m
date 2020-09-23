clear all;
close all
%example plots of different umap results

%plot unsupervised
load([pwd '\umap_data\umap_unsupervised_cluster_perc.mat'])
figure;
for pl=1:size(perc,2)
    subplot(3,3,pl);bar(perc{pl}')
    set(gca,'xlim',[-.5 (size(perc{pl},2))+.5])
    if pl==4
        ylabel('Number of snippets (%)')
    end
    if pl==7 || pl==8
        xlabel('Cluster #')
    end
    title(['Bin: ' num2str(bin(pl))])
    if pl==1
        legend('Saline','Deaf')
        legend boxoff;
    end
end
load([pwd '\umap_data\pdfdiffs.mat'])
figure;
plot(bin,dval,'*k')
set(gca,'xlim',[bin(1)-5 bin(end)+5],'xtick',bin,'ylim',[min(dval)-.1*max(dval) max(dval)+.1*max(dval)])
xlabel('Bin (ms)')
ylabel('PDF difference (abs(sum(pdfsal-pdfdeaf))')


%plot supervised
load([pwd '\umap_data\umap_supervised_perc.mat'])
figure;
for pl=1:size(dvals,1)
    pcs=squeeze(dvals(pl,:,:));
    [nbins,indx]=hist(pcs(:),0:10:100);
    subplot(3,3,pl);bar(indx,nbins./sum(nbins).*100)
    set(gca,'ylim',[0 50],'xlim',[indx(1)-.5 indx(end)+.5],'xtick',(0:25:100))
    if pl==4
        ylabel('Number classification (%)')
    end
    if pl==7 || pl==8
        xlabel('Correct classification (%)')
    end
    title(['Bin: ' num2str(bin(pl))])
end

