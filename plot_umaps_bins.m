clear all;
close all;
%examples of the unsupervised umap projections for different bins
load([pwd '\umap_data\umap_data_saldeaf_spec_tsc.mat'])
for bnr=1:length(bin)
    datadeaf=prepDATAd.umap{bnr};
    datasal=prepDATAs.umap{bnr};
    iter=size(datasal,1);
    MAT=[datasal;datadeaf];
    indx=[];
    indx{1}=1:iter;
    indx{2}=iter+1:size(MAT,1);
    clear MAT datadeaf datasal;
    load([pwd '\umap_data\UMAP_saldeaf_' num2str(bin(bnr)) '_tsc.mat'])
    figure(1);subplot((length(bin)/2),length(bin)/2,bnr+floor((bnr-1)/4)*4)
    plot(reduction(indx{1},1),reduction(indx{1},2),'k.','markersize',5)
    title(['Bin ' num2str(bin(bnr)) ', Saline'])
    figure(1);subplot((length(bin)/2),length(bin)/2,bnr+floor((bnr-1)/4)*4+4)
    plot(reduction(indx{2},1),reduction(indx{2},2),'r.','markersize',5)
        title(['Bin ' num2str(bin(bnr)) ', Deaf'])
end
