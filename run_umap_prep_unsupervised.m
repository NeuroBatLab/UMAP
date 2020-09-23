clear all;
close all;
%example of unsupervised umap on time-scaled data prepped through
%prep_audio_umap.m; runs unsupervised UMAP, no template for each existing
%bin and in this case for a saline and deaf data set with generated ID-set

%set plotting of results (=1)
plotdata=0;
%run umap(=1) or grab previous run (=0)
run=1;
%calculate the probablity density function between umap projections (=1)
runpdf=0;
%example of prepped data
load([pwd '\umap_data\umap_data_saldeaf_spec_tsc.mat'])
for bnr=1:length(bin)
    disp(['Running bin ' num2str(bin(bnr)) '...'])
    datadeaf=prepDATAd.umap{bnr};
    datasal=prepDATAs.umap{bnr};
    iter=size(datasal,1);
    MAT=[datasal;datadeaf];
    indx=[];
    indx{1}=1:iter;
    indx{2}=iter+1:size(MAT,1);
    biossal=prepDATAs.bios{bnr};
    biosdeaf=prepDATAd.bios{bnr};
    %generates ID for each snippet
    for ids=1:2
        if ids==1
            Pos=prepDATAs.Pos{bnr};
            IDs=prepDATAs.IDs;
        else
            Pos=prepDATAd.Pos{bnr};
            IDs=prepDATAd.IDs;
        end
        newIDs=[];
        for i=1:length(Pos)
            newIDs=[newIDs ones(1,Pos(i)).*IDs(i)];
        end
        allIDs{ids}=newIDs;
    end
    if run
        addpath([pwd '\umap\'])
        addpath([pwd '\umap\sgdCpp_files'])
        addpath([pwd '\util\'])
        [reduction, umap, clustering]=run_umap(MAT,'verbose','none');
        save([pwd '\umap_data\UMAP_saldeaf_' num2str(bin(bnr)) '_tsc.mat'],...
            'umap','reduction','clustering','indx','-v7.3');
    else
        load([pwd '\umap_data\UMAP_saldeaf_' num2str(bin(bnr)) '_tsc.mat']);
    end
    %compare clustering between saline and deaf
    pc=[];
    for i=1:2
        [nbins,n]=hist(clustering(indx{i}),0:1:double(max(clustering)));
        pc(i,:)=nbins./length(indx{i}).*100;
    end
    perc{bnr}=pc;
    %runs and plots the probablity density functions between saline and
    %deaf projections
    if runpdf
        N=20;
        [rx,ry,rpdfs]=jointpdf(reduction(indx{1},1),reduction(indx{1},2),N,plotdata);
        [rx,ry,rpdfd]=jointpdf(reduction(indx{2},1),reduction(indx{2},2),N,plotdata);
        pdfdiff=(rpdfs)-(rpdfd);
        if plotdata==1
            figure;surf(rx,ry,pdfdiff)
        end
        %absolute difference value of pdfs
        dval(bnr)=abs(sum(sum(pdfdiff)));
    end
    if plotdata==1
        %plots projection of saline and deaf separately
        figure;colorstr={'k';'r'};tits={'Saline';'Deaf'};
        for pl=1:2
            subplot(1,2,pl)
            plot(reduction(indx{pl},1),reduction(indx{pl},2),[colorstr{pl} '.'],'Markersize',5)
            
            title(['Bin' num2str(bin(bnr)) ': ' tits{pl}])
        end
        xlims=xlim;
        ylims=ylim;
        set(gca,'ylim',[ylims],'xlim',[xlims])
        %color codes the saline or deaf projection according to the biosound
        %parameters
        titstr={'tAmp';'Amp';'Sal';'F0';'SpecM'};
        for bpl=1:5
            BIO=biossal(:,bpl);
            BIO=[BIO;biosdeaf(:,bpl)];
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
                title(['Bin' num2str(bin(bnr)) ': ' tits{pl} '-' titstr{bpl}])
                set(gca,'ylim',[ylims],'xlim',[xlims])
            end
        end
        %plots individual saline bat or deaf projections
        for pl=1:2
            xData=reduction(indx{pl},1);
            yData=reduction(indx{pl},2);
            batIDs=allIDs{pl};
            batid=unique(batIDs);
            figure;
            for bn=1:length(batid)
                indxid=find(batIDs==batid(bn));
                subplot(3,2,bn)
                hScat = scatter(xData(indxid), yData(indxid), 20,colorstr{pl},'Filled');
                set(gca,'ylim',[ylims],'xlim',[xlims])
                title(['Bin' num2str(bin(bnr)) ': ' tits{pl} '-Ind.' num2str(bn)])
            end
        end
        %plots individual saline or deaf bat projections, color coded by
        %the biosound parameters
        for pl=1:2
            batIDs=allIDs{pl};
            batid=unique(batIDs);
            for bpl=1:5
                figure;
                for bn=1:length(batid)
                    xData=reduction(indx{pl},1);
                    yData=reduction(indx{pl},2);
                    indxid=find(batIDs==batid(bn));
                    xData=xData(indxid);
                    yData=yData(indxid);
                    if pl==1
                        BIO=biossal(indxid,bpl);
                    else
                        BIO=biosdeaf(indxid,bpl);
                    end
                    subplot(3,2,bn)
                    hScat = scatter(xData, yData, 20, BIO, 'Filled');
                    hcbar = colorbar;
                    hcbar.Label.String='Bio value';
                    caxis([allax(bpl,1) allax(bpl,2)])
                    title(['Bin' num2str(bin(bnr)) ': ' tits{pl} '-' titstr{bpl} '-B' num2str(bn)])
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
            title(['Bin' num2str(bin(bnr)) ': ' tits{pl}])
        end
    end
end
if runpdf
    save([pwd '\umap_data\pdfdiffs.mat'],'dval','bin','N')
end
%save([pwd '\umap_data\umap_unsupervised_cluster_perc.mat'],'perc','bin')
