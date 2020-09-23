clear all;
close all;
%example of supervised umap on time-scaled data prepped through
%prep_audio_umap.m; runs supervised UMAP, template (training data-set) is generated from defined
%number of vocalizations for each saline individual for each existing
%bin and for  fixed min_dis and n_neighbor values retrieved through intial
%clustering testing

addpath([pwd '\umap\'])
addpath([pwd '\umap\sgdCpp_files'])
addpath([pwd '\util\'])
%best min_dis and n_neighbor values from intial run
min_dis=0.051;
n_neighbors =3;
%minimal data for umap
minN_umap=100;
%maximal valid projection range for umap
maxproj_umap=500;
%to grab all data for calculation set to 0, else 1:
randomize_data=1;
%proportion of data to grab randomly if independent of individual
randprop=.4;
%randomize weighted among individuals (=same number of calls per
%individual) set to 1, else set 0
rand_weighted=1;
%# of repetitions for randomly grabbing N calls for each individual for data-sets
rep=10;
%plot results (=1)
plotdata=0;
%example of prepped data
load([pwd '\umap_data\Julie_umap_data_saldeaf_spec_tsc.mat'])
cluscompd=NaN(length(bin),rep,2,100);
cluscomph=NaN(length(bin),rep,2,100);
%run supervised UMAP
for bnr=1:length(bin)
    disp(['Bin ' num2str(bin(bnr))])
    datadeaf=prepDATAd.umap{bnr};
    datasal=prepDATAs.umap{bnr};
    biosdeaf=prepDATAd.bios{bnr};
    biossal=prepDATAs.bios{bnr};
    %generates ID for each snippet
    allindx=NaN(2,5);
    for ids=1:2
        if ids==1
            IDs=prepDATAs.IDs;
        else
            IDs=prepDATAd.IDs;
        end
        batid=unique(IDs);
        if ids==1
            indxidsall=batid;
        else
            indxiddall=batid;
        end
        for bn=1:length(batid)
            indxid(bn)=length(find(IDs==batid(bn)));
        end
        allindx(ids,1:length(indxid))=indxid;
    end
    %define maximal possible number of calls for each individual
    callN=min(min(allindx));
    %run through the different train/test combinations
    if randomize_data
        %grab the random calls (per individual if rand_weighted = 1)
        if rand_weighted~=1
            IDs=prepDATAs.IDs;
            callN=round(randprop*length(IDs));
        end
        [sDATAtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,indxidsall,callN,rep,rand_weighted);
        if rand_weighted~=1
            IDs=prepDATAd.IDs;
            callN=round(randprop*length(IDs));
        end
        [dDATAtrindx]=get_random_calls(prepDATAd.Pos{bnr},prepDATAd.IDs,indxiddall,callN,rep,rand_weighted);
    else
        %grab all data
        rep=1;
        sDATAtrindx{1}=1:size(datasal,1);
        dDATAtrindx{1}=1:size(datadeaf,1);
        cluscompd=NaN(length(bin),rep,2,100);
        cluscomph=NaN(length(bin),rep,2,100);
    end
    %run the number repetitions
    for rp=1:rep
        %if not enough data for umap
        if length(sDATAtrindx{rp})<minN_umap
            disp('Nsal too small for this rep! Redoing random data set generation independent of individual')
            IDs=prepDATAs.IDs;
            callN=round(randprop*length(IDs));
            [newtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,indxidsall,callN,1,0);
            sDATAtrindx{rp}=newtrindx{1};
            %if not enough data for umap
            if length(sDATAtrindx{rp})<minN_umap
                disp('Nsal too small! Grabbing more data')
                callN=round(.8*length(IDs));
                [newtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,indxidsall,callN,1,0);
                sDATAtrindx{rp}=newtrindx{1};
                %if not enough data for umap
                if length(sDATAtrindx{rp})<minN_umap
                    disp('Nsal too small! Grabbing all data')
                    sDATAtrindx{rp}=1:size(datasal,1);
                end
            end
        end
        if length(dDATAtrindx{rp})<minN_umap
            disp('Ndeaf too small for this rep! Redoing random data set generation independent of individual')
            IDs=prepDATAd.IDs;
            callN=round(randprop*length(IDs));
            [newtrindx]=get_random_calls(prepDATAd.Pos{bnr},prepDATAd.IDs,indxiddall,callN,1,0);
            dDATAtrindx{rp}=newtrindx{1};
            %if not enough data for umap
            if length(dDATAtrindx{rp})<minN_umap
                disp('Ndeaf too small! Grabbing more data')
                callN=round(.8*length(IDs));
                [newtrindx]=get_random_calls(prepDATAd.Pos{bnr},prepDATAd.IDs,indxiddall,callN,1,0);
                dDATAtrindx{rp}=newtrindx{1};
                %if not enough data for umap
                if length(dDATAtrindx{rp})<minN_umap
                    disp('Ndeaf too small! Grabbing all data')
                    dDATAtrindx{rp}=1:size(datadeaf,1);
                end
            end
        end
        disp(['Rep #' num2str(rp) ' of ' num2str(rep)])
        Train=datasal(sDATAtrindx{rp},:);
        Test=datadeaf(dDATAtrindx{rp},:);
        %generate labels
        if rp==1
            labs={size(Train,2)};
            for f=1:size(Train,2)
                labs{f}=['f' num2str(f)];
            end
        end
        disp(['Running umaps...'])
        %run template=train data-set
        [reduction, umap,clustering]=run_umap(Train,'parameter_names',labs,...
            'min_dis',min_dis,'n_neighbors',n_neighbors,...
            'verbose','none');
        figure(8000);plot(reduction(:,1),reduction(:,2),'.','Markersize',5)
        xlims=xlim;
        %if umap projection is ok
        if max(abs(xlims))<maxproj_umap
            save([pwd '\umap_data\UMAP_saltrain_tsc.mat'], 'umap','reduction','-v7.3');
            %run test with template=train data-set
            [reductiont, umapt,clusteringt]= run_umap(Test, 'template_file',...
                [pwd '\umap_data\UMAP_saltrain_tsc.mat'],'parameter_names',labs,...
                'verbose','none');
            disp(['Running clustering...'])
            %dbscan
            vec=[reduction;reductiont];
            clindx{1}=1:size(reduction,1);
            clindx{2}=size(reduction,1)+1:size(vec,1);
            D=pdist2(vec,vec);
            [idx, corepts] = dbscan(D,1,5,'Distance','precomputed');
            tits={'training';'test'};
            
            for mi=1:2
                if max(idx)>0
                    [nbins,idxh]=hist(idx(clindx{mi}),-1:1:max(idx));
                    cluscompd(bnr,rp,mi,1:length(nbins))=nbins./length(clindx{mi});
                end
                if plotdata
                    
                    if mi==1
                        red=reduction;
                    else
                        red=reductiont;
                    end
                    figure(bnr);subplot(2,1,mi);
                    plot(red(:,1),red(:,2),'.','Markersize',5)
                    title(['Bin :' num2str(bin(bnr))])
                    if mi==1
                        xlims=xlim;
                        ylims=ylim;
                    end
                    figure(100+bnr);subplot(2,2,mi)
                    hScat = scatter(vec(clindx{mi},1), vec(clindx{mi},2), 20, idx(clindx{mi}), 'Filled');
                    title(['Bin: ' num2str(bin(bnr)) '-dbscan: ' tits{mi}])
                    set(gca,'xlim',xlims,'ylim',ylims)
                end
            end
            
            %hdbscan
            clusteringh=HDBSCAN(vec);
            clusteringh.minpts=10;
            clusteringh.minclustsize=4;
            fit_model(clusteringh)
            get_best_clusters(clusteringh)
            get_membership(clusteringh)
            idx=double(clusteringh.labels);
            
            for mi=1:2
                if max(idx)>0
                    [nbins,idxh]=hist(idx(clindx{mi}),-1:1:max(idx));
                    cluscomph(bnr,rp,mi,1:length(nbins))=nbins./length(clindx{mi});
                end
                if plotdata
                    figure(100+bnr);subplot(2,2,mi+2)
                    hScat = scatter(vec(clindx{mi},1), vec(clindx{mi},2), 20, idx(clindx{mi}), 'Filled');
                    title(['Bin: ' num2str(bin(bnr)) '-hbscan: ' tits{mi}])
                    set(gca,'xlim',xlims,'ylim',ylims)
                    
                    figure(2000+bnr+100*(mi-1));
                    %color codes the saline or deaf projection according to the biosound
                    %parameters
                    titstr={'tAmp';'Amp';'Sal';'F0';'SpecM'};
                    for bpl=1:5
                        if mi==1
                            BIO=biossal(:,bpl);
                        else
                            BIO=biosdeaf(:,bpl);
                        end
                        subplot(3,2,bpl);
                        hScat = scatter(vec(clindx{mi},1), vec(clindx{mi},2), 20, BIO, 'Filled');
                        hcbar = colorbar;
                        hcbar.Label.String='Bio value';
                        caxis([min(BIO) max(BIO)])
                        allax(bpl,:)=[min(BIO) max(BIO)];
                        title(['Bin' num2str(bin(bnr)) ': ' tits{mi} '-' titstr{bpl}])
                        set(gca,'ylim',[ylims],'xlim',[xlims])
                        
                    end
                end
            end
            
        else
            disp('umap projection not successful...')
        end
        close(figure(8000));
    end
    
    
end

count=0;
rp=5;
for bnr=1:length(bin)
    alli=[];
    for mi=1:2
        indx=find(~isnan(cluscompd(bnr,rp,mi,:)));
        if ~isempty(indx)
            alli(mi)=indx(end);
        end
    end
    if ~isempty(alli)
        count=count+1;
        figure(2000);subplot(3,3,count);
        bar(1:max(alli),squeeze(cluscompd(bnr,rp,:,1:max(alli)))')
        title(['Bin: ' num2str(bin(bnr))])
        if count==1
            legend('Control','Deaf')
        end
        if bnr==1
            
        elseif count==4
            ylabel('Cluster (%)')
        elseif count==5
            xlabel('Cluster #')
        end
    end
end