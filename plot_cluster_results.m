function plot_cluster_results(DATA,BIOS,idx,bio_indx,min_dis,n_neighbors,plotpos,bin,LIMS,bioflag,fh1,fh2)
%plot dbscan and hdbscan clustering results
figure(fh1);subplot(ceil(length(n_neighbors)/ceil(length(n_neighbors)/2)),...
    ceil(length(n_neighbors)/2),plotpos.nn);
hScat = scatter(DATA(:,1), DATA(:,2), 20, idx.d, 'Filled');
title(['NN: ' num2str(n_neighbors(plotpos.nn))])
figure(fh2);subplot(ceil(length(n_neighbors)/ceil(length(n_neighbors)/2)),...
    ceil(length(n_neighbors)/2),plotpos.nn);
hScat = scatter(DATA(:,1), DATA(:,2), 20, idx.h, 'Filled');
title(['NN: ' num2str(n_neighbors(plotpos.nn))])
if bioflag
    %color codes the projection according to the biosound
    %parameters
    titstr={'tAmp';'Amp';'Sal';'F0';'SpecM'};
    figure('Name',['Control_init-Biosound: Bin' num2str(bin(plotpos.bnr)) ': MD' num2str(min_dis(plotpos.md))...
            '_NN' num2str(n_neighbors(plotpos.nn)) '-Rep#' num2str(plotpos.rep)],'NumberTitle','off');
    for bpl=1:5
        subplot(3,2,bpl);
        hScat = scatter(DATA(:,1), DATA(:,2), 20, BIOS(bio_indx,bpl), 'Filled');
        hcbar = colorbar;
        hcbar.Label.String='Bio value';
        caxis([min(BIOS(:,bpl)) max(BIOS(:,bpl))])
        title([titstr{bpl}])
        set(gca,'ylim',[LIMS.ylims],'xlim',[LIMS.xlims])
    end
end