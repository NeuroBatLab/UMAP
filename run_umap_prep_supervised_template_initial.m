clear all;
close all;
%example of supervised umap on time-scaled data prepped through
%prep_audio_umap.m; runs supervised UMAP, template (training data-set)is generated from defined
%number of vocalizations for each saline individual for each existing
%bin and for different min_dis and n_neighbor values
%clustering quality is evaluated through Hopkins method

addpath([pwd '\umap\'])
addpath([pwd '\umap\sgdCpp_files'])
addpath([pwd '\util\'])
%different min_dis and n_neighbor values
min_dis=[0.051 .1 .2];% .4 .6 .79];
n_neighbors =[3 5 10];% 20];% 50 100 150 199];
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
%# of repetitions for randomly grabbing N calls for each individual for
%training data-set
rep=10;
%distribution of random sample data for hopkins evaluation
distribution='uniform_convex_hull';
%neighbor value for hopkins stats
n_neighbor=1;
%proportion of umap data into hopkins
prop=.8;
%proportion of umap data not used for hopkins evaluation
portion=0.05;
%# of repetitions for hopkins calculation on percentage (=prop) of umap data
%-prop of that data
rephop=20;
%plot results (=1)
plotdata=0;
%example of prepped data
load([pwd '\umap_data\Julie_umap_data_saldeaf_spec_tsc.mat'])
%run supervised UMAP
for bnr=1:length(bin)
    disp(['Bin ' num2str(bin(bnr))])
    datasal=prepDATAs.umap{bnr};
    biossal=prepDATAs.bios{bnr};
    %generates ID for each snippet
    IDs=prepDATAs.IDs;
    batid=unique(IDs);
    allindx=zeros(1,length(batid));
    for bn=1:length(batid)
        allindx(bn)=length(find(IDs==batid(bn)));
    end
    %define maximal possible number of calls for each individual
    callN=min(min(allindx));
    if rand_weighted~=1
        callN=round(randprop*length(IDs));
    end
    if randomize_data
        %grab the random calls (per individual if rand_weighted = 1)
        [sDATAtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,batid,callN,rep,rand_weighted);
    else
        %grab all data
        rep=1;
        sDATAtrindx{1}=1:size(datasal,1);
    end
    if bnr==1
        clN=NaN(length(bin),rep,length(min_dis),length(n_neighbors),2);
        hops=NaN(length(bin),rep,length(min_dis),length(n_neighbors),3);
    end
    for rp=1:rep
        %if not enough data for umap
        if length(sDATAtrindx{rp})<minN_umap
            disp('N too small for this rep! Redoing random data set generation independent of individual')
            callN=round(randprop*length(IDs));
            [newtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,batid,callN,1,0);
            sDATAtrindx{rp}=newtrindx{1};
            %if not enough data for umap
            if length(sDATAtrindx{rp})<minN_umap
                disp('N too small! Grabbing more data')
                callN=round(.8*length(IDs));
                [newtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,batid,callN,1,0);
                sDATAtrindx{rp}=newtrindx{1};
                %if not enough data for umap
                if length(sDATAtrindx{rp})<minN_umap
                    disp('N too small! Grabbing all data')
                    sDATAtrindx{rp}=1:size(datasal,1);
                end
            end
        end
        disp(['Rep #' num2str(rp) ' of ' num2str(rep)])
        Train=datasal(sDATAtrindx{rp},:);
        %generate labels
        if rp==1
            labs={size(Train,2)};
            for f=1:size(Train,2)
                labs{f}=['f' num2str(f)];
            end
        end
        for md=1:length(min_dis)
            for nn=1:length(n_neighbors)
                disp(['MD: ' num2str(min_dis(md)) ' - NN: ' num2str(n_neighbors(nn))])
                disp(['Running umap...'])
                %run template=train data-set
                [reduction, umap,clustering]=run_umap(Train,'parameter_names',labs,...
                    'min_dis',min_dis(md),'n_neighbors',n_neighbors(nn),...
                    'verbose','none');
                figure(8000);plot(reduction(:,1),reduction(:,2),'.','Markersize',5)
                drawnow;
                xlims=xlim;
                %if umap projection is ok
                if max(abs(xlims))<maxproj_umap
                    disp(['Running clustering...'])
                    %dbscan
                    D=pdist2(reduction,reduction);
                    [idx, corepts] = dbscan(D,1,5,'Distance','precomputed');
                    clN(bnr,rp,md,nn,1)=length(unique(idx));
                    if plotdata
                        figure(md);subplot(2,2,nn);
                        hScat = scatter(reduction(:,1), reduction(:,2), 20, idx, 'Filled');
                        title(['dbscan: MD: ' num2str(min_dis(md)) '_NN: ' num2str(n_neighbors(nn))])
                    end
                    %hdbscan
                    clusteringh=HDBSCAN(reduction);
                    clusteringh.minpts=5;
                    clusteringh.minclustsize=5;
                    fit_model(clusteringh)
                    get_best_clusters(clusteringh)
                    get_membership(clusteringh)
                    idx=clusteringh.labels;
                    clN(bnr,rp,md,nn,2)=length(unique(idx));
                    if plotdata
                        figure(md+100);subplot(2,2,nn);
                        hScat = scatter(reduction(:,1), reduction(:,2), 20, idx, 'Filled');
                        title(['hdbscan: MD: ' num2str(min_dis(md)) '_NN: ' num2str(n_neighbors(nn))])
                    end
                    disp(['Running hopkins...'])
                    %hopkins cluster evaluation
                    H=NaN(rephop,1);
                    for rh=1:rephop
                        seq=randperm(size(reduction,1));
                        seq=seq(1:round(size(reduction,1)*prop));
                        H(rh)=hopkins_stats(reduction(seq,:),min(reduction(seq,:),[],1), ...
                            max(reduction(seq,:),[],1),round(prop*size(reduction,1))-...
                            round(portion*(prop*size(reduction,1))),distribution,n_neighbor,flip);
                    end
                    %bootstrap hopkins results and get CI
                    CI=bootci(1e3,@mean,H);
                    hops(bnr,rp,md,nn,:)=[mean(H);CI];
                else
                    disp('umap projection not successful...')
                end
                close(figure(8000));
            end
        end
    end
end
save([pwd '\umap_data\Julie_umap_data_saldeaf_spec_tsc_salineumapinit.mat'],'hops','clN')
for bnr=1:size(hops,1)
    count=0;
    for md=1:size(hops,3)
        for nn=1:size(hops,4)
            count=count+1;
            figure(bnr);hold on;
            plot(count,squeeze(hops(bnr,:,md,nn,1)),'k*')
            %plot(count,squeeze(hops(bnr,:,md,nn,2)),'r*')
            %plot(count,squeeze(hops(bnr,:,md,nn,3)),'g*')
            xlab{count}=[num2str(md) '_' num2str(nn)];
        end
    end
    set(gca,'xlim',[0 count+1],'xtick',(1:1:count),'xticklabel',xlab)
end
