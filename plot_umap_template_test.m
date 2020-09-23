function plot_umap_template_test(DATA,VEC,IDX,plindx,bin,minmxidx,corval,BIOS,bios,fhigh,plotpos,LIMS,bioflag,fh)


tits={'training';'test'};
%plot clustering of umap reduction
figure(fh);subplot(2,1,plotpos.mi)
hScat = scatter(DATA(plindx,1), DATA(plindx,2), 20,IDX(plindx), 'Filled');
caxis([min(IDX) max(IDX)])
hcbar = colorbar;
hcbar.Label.String='Cluster #';
if plotpos.mi==1
title([tits{plotpos.mi}])
else
 title([tits{plotpos.mi} '-Comp#' num2str(plotpos.compN)])   
end
set(gca,'xlim',LIMS.xlims,'ylim',LIMS.ylims)
%revert and plot spectrograms of each cluster
for cln=minmxidx(1): minmxidx(2)
    clnindx=find(IDX(plindx)==cln);
    plot_spec_cluster(VEC(plindx(clnindx)-corval,:),...
        bin(plotpos.bnr),fhigh,plotpos.compN,cln,plotpos.mi,plotpos.rep)
end
if bioflag
    if plotpos.mi==1
        figure('Name',['Biosound: Bin' num2str(bin(plotpos.bnr)) '-' tits{plotpos.mi} '-Rep#' num2str(plotpos.rep)],'NumberTitle','off');
    else
        figure('Name',['Biosound: Bin' num2str(bin(plotpos.bnr)) '-' tits{plotpos.mi} '-Comp#' num2str(plotpos.compN) '-Rep#' num2str(plotpos.rep)],'NumberTitle','off');
    end
    %color codes the saline or deaf projection according to the biosound
    %parameters
    titstr={'tAmp';'Amp';'Sal';'F0';'SpecM'};
    for bpl=1:5
        BIO=BIOS(:,bpl);
        biovec=bios(:,bpl);
        subplot(3,2,bpl);
        hScat = scatter(DATA(plindx,1), DATA(plindx,2), 20, biovec, 'Filled');
        hcbar = colorbar;
        hcbar.Label.String='Bio value';
        caxis([min(BIO) max(BIO)])
        title([titstr{bpl}])
        set(gca,'ylim',[LIMS.ylims],'xlim',[LIMS.xlims])
    end
end