clear all;
close all;
%example of unsupervised umap on dynamic time-warped data prepped through
%prep_audio_umap.m; runs unsupervised UMAP, no template

%set plotting of results (=1)
plotdata=1;
%run umap(=1) or grab previous run (=0)
run=0;
%saline
IDs{1}=[11648;14461;14464;65696;71353];
%deaf
IDs{2}=[14463;71043;71047;71351;71354];
%example prepped data
load([pwd '\umap_data\umap_data_saldeaf_spec_dtw.mat'])
if run
    addpath([pwd '\umap\'])
    addpath([pwd '\umap\sgdCpp_files'])
    addpath([pwd '\util\'])
    min_dis=[0.051 .1 .4 .6 .79];
    n_neighbors =[3 5 10 20];
    for loop1=1:length(min_dis)
        figure;
        for loop2=1:length(n_neighbors)
            [reduction, umap,clustering]=run_umap(prepDATA.umap{1},'metric','precomputed',...
                'verbose','none','n_neighbors',n_neighbors(loop2),'min_dis',min_dis(loop1));
            for pl=1:2
                subplot(2,4,loop2+(pl-1)*4);plot(reduction(indx{pl},1),reduction(indx{pl},2),...
                    '.','color',colorstr{pl},'Markersize',5)
                title([tits{pl} '_MD:' num2str(min_dis(loop1)) '_NN:' num2str(n_neighbors(loop2))])
            end
        end
    end
    %     save([pwd '\umap_data\UMAP_saldeaf_dtw.mat'],...
    %         'umap','reduction','clustering','-v7.3');
else
    %load([pwd '\umap_data\UMAP_saldeaf_dtw.mat'])
end
for pl=1:2
    ind=[];
    ID=IDs{pl};
    for i=1:length(ID)
        ind=[ind;find(prepDATA.IDs==ID(i))];
    end
    indx{pl}=ind;
end
if run
    MAT=prepDATA.umap{1};
    Train=MAT([1:ceil(length(indx{1})/2) length(indx{1})+1:floor(length(indx{2})/2)+length(indx{1})],...
        [1:ceil(length(indx{1})/2) length(indx{1})+1:floor(length(indx{2})/2)+length(indx{1})]);
    Test=MAT([ceil(length(indx{1})/2)+1:length(indx{1}) length(indx{1})+floor(length(indx{2})/2)+1:end],...
        [ceil(length(indx{1})/2)+1:length(indx{1}) length(indx{1})+floor(length(indx{2})/2)+1:end]);
    Train(:,end+1)=0;
    Train(ceil(length(indx{1})/2)+1:end,end)=1;
    labs=[];
    for f=1:size(Train,2)-1
        l=['f' num2str(f)];
        labs{f}=l;
    end
    %run template=train data-set
    [reduction, umap,clustering]=run_umap(Train,'parameter_names',labs,...
        'label_column',size(Train,2),'metric','precomputed',...
        'verbose','graphic');
    save([pwd '\umap_data\UMAP_saldeaftrain_dtw.mat'], 'umap','reduction','clustering','-v7.3');
    
    [reductiont, umapt,clustering]= run_umap(Test, 'template_file',...
        [pwd '\umap_data\UMAP_saldeaftrain_dtw.mat'],'parameter_names',labs,...
        'metric','precomputed','verbose','graphic');
end
if plotdata
    %plots projection of saline and deaf separately
    colorstr={'k';'r'};tits={'Saline';'Deaf'};
    figure;
    for pl=1:2
        subplot(2,1,pl);plot(reduction(indx{pl},1),reduction(indx{pl},2),...
            '.','color',colorstr{pl},'Markersize',5)
        title(tits{pl})
    end
    xlims=xlim;
    ylims=ylim;
    %color codes the saline or deaf projection according to the mean biosound
    %parameters of each vocalization
    allbios=prepDATA.bios{1};
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
    %plots individual saline bat or deaf projections
    for pl=1:2
        xData=reduction(indx{pl},1);
        yData=reduction(indx{pl},2);
        batIDs=prepDATA.IDs(indx{pl});
        batid=unique(batIDs);
        figure;
        for bn=1:length(batid)
            indxid=find(batIDs==batid(bn));
            subplot(3,2,bn)
            hScat = scatter(xData(indxid), yData(indxid), 20,colorstr{pl},'Filled');
            set(gca,'ylim',[ylims],'xlim',[xlims])
            title([tits{pl} '-Ind.' num2str(bn)])
        end
    end
    %plots individual saline or deaf bat projections, color coded by
    %the biosound parameters
    for pl=1:2
        batIDs=prepDATA.IDs(indx{pl});
        batid=unique(batIDs);
        bios=allbios(indx{pl},:);
        for bpl=1:5
            figure;
            for bn=1:length(batid)
                xData=reduction(indx{pl},1);
                yData=reduction(indx{pl},2);
                indxid=find(batIDs==batid(bn));
                xData=xData(indxid);
                yData=yData(indxid);
                BIO=bios(indxid,bpl);
                subplot(3,2,bn)
                hScat = scatter(xData, yData, 20, BIO, 'Filled');
                hcbar = colorbar;
                hcbar.Label.String='Bio value';
                caxis([allax(bpl,1) allax(bpl,2)])
                title([tits{pl} '-' titstr{bpl} '-B' num2str(bn)])
                set(gca,'ylim',[ylims],'xlim',[xlims])
            end
        end
    end
    %color codes the saline or deaf projections based on the snippet an
    %call sequence
    figure;
    for pl=1:2
        xData=reduction(indx{pl},1);
        yData=reduction(indx{pl},2);
        cVect=1:length(xData);
        subplot(1,2,pl)
        hScat = scatter(xData, yData, 20, cVect, 'Filled');
        hcbar = colorbar;
        hcbar.Label.String='Row #';
        caxis([1 length(xData)])
        title([tits{pl}])
    end
end