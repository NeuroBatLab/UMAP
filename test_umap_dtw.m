
clear all;
close all;
pn='C:\Users\Daria\audio_stuff\umap_data\upmap_data_controldeaf_spec_ed.mat';
load(pn)
neigh=10;
min_d=.7;

[reduction, umap,clustering]=run_umap(EDs, 'n_neighbors',neigh,...
    'min_dist',min_d,'metric','precomputed');
xlims=xlim;
ylims=ylim;
for run=1:2
    if run==2
        %deaf
        IDs=[14463;71043;71047;71351;71354];
    else
        %saline
        IDs=[11648;14461;14464;65696;71353];
    end
    ind=[];
    for i=1:length(IDs)
        ind=[ind;find(batIDs==IDs(i))];
    end
    indx{run}=ind;
end
colorstr={'k';'r'};tits={'Saline';'Deaf'};
figure;
for pl=1:2
subploy(2,1,pl);plot(reduction(indx{pl},1),reduction(indx{pl},2),'.','color',colorstr{pl},'Markersize',5)
title(tits{pl})
end
load(['C:\Users\Daria\audio_stuff\umap_data\upmap_data_controldeaf_bio_v2_dtw.mat'])
titstr={'tAmp';'Amp';'Sal';'F0';'SpecM'};
for bpl=1:5
    BIO=allbios(:,bpl);
    figure;
    for pl=1:2
        xData=reduction(indx{pl},1);
        yData=reduction(indx{pl},2);
        subplot(1,2,pl)
        hScat = scatter(xData, yData, 20, BIO(indx{pl}), 'Filled');
        hcbar = colorbar;
        hcbar.Label.String='Bio value';
        caxis([min(BIO) max(BIO)])
        allax(bpl,:)=[min(BIO) max(BIO)];
        title([tits{pl} '-' titstr{bpl}])
        set(gca,'ylim',[ylims],'xlim',[xlims])
    end
end
