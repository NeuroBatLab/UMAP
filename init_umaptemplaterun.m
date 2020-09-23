function [hops,clN]=init_umaptemplaterun(prepDATA,bin,PARAS,plotdata,bioflag)
%initial run of umap on time-scaled data prepped through
%prep_audio_umap.m to identify optimal minimal distance and nearest neighbor pair and
%ideal clustering method;
%runs supervised UMAP, template (training data-set)is generated from defined
%number of vocalizations for each saline individual for each existing
%bin and for different min_dis and n_neighbor values
%clustering quality is evaluated through Hopkins method
%requires umap path to be set and
%functions needed: get_random_calls, check_dataforumap, run_clustering, dbscan,
%HDBSCAN, hopkins_stats

%INPUT: prepDATA, PARAS, plotdata(0 or 1)
%OUTPUT: hopkin evaluations and cluster numbers

%different min_dis and n_neighbor values
min_dis=PARAS.min_dis;
n_neighbors =PARAS.n_neighbors;
%minimal data for umap
minN_umap=PARAS.minN_umap;
%maximal valid projection range for umap
maxproj_umap=PARAS.maxproj_umap;
%to grab all data for calculation set to 0, else 1:
randomize_data=PARAS.randomize_data;
%proportion of data to grab randomly if independent of individual
randprop=PARAS.randprop;
%randomize weighted among individuals (=same number of calls per
%individual) set to 1, else set 0
rand_weighted=PARAS.rand_weighted;
%# of repetitions for randomly grabbing N calls for each individual for
%training data-set
rep=PARAS.rep;
%distribution of random sample data for hopkins evaluation
distribution=PARAS.distribution;
%flip hopkins
flip=PARAS.flip;
%neighbor value for hopkins stats
n_neighborhop=PARAS.n_neighborhop;
%proportion of umap data into hopkins
prop=PARAS.prop;
%proportion of umap data not used for hopkins evaluation
portion=PARAS.portion;
%# of repetitions for hopkins calculation on percentage (=prop) of umap data
%-prop of that data
rephop=PARAS.rephop;
db_paras=PARAS.db_paras;
hdb_paras=PARAS.hdb_paras;
%run supervised UMAP
for bnr=1:length(bin)
    disp(['Bin ' num2str(bin(bnr))])
    data=prepDATA.umap{bnr};
    bios=prepDATA.bios{bnr};
    %generates ID for each snippet
    IDs=prepDATA.IDs;
    if ~isempty(IDs)
        batid=unique(IDs);
        allindx=zeros(1,length(batid));
        for bn=1:length(batid)
            allindx(bn)=length(find(IDs==batid(bn)));
        end
        %define maximal possible number of calls for each individual
        callN=min(min(allindx));
    end
    if rand_weighted~=1
        callN=round(randprop*length(IDs));
    end
    if randomize_data
        %grab the random calls (per individual if rand_weighted = 1)
        [DATAtrindx]=get_random_calls(prepDATA.Pos{bnr},IDs,batid,callN,rep,rand_weighted);
    else
        %grab all data
        rep=1;
        DATAtrindx{1}=1:size(data,1);
    end
    if bnr==1
        clN=NaN(length(bin),rep,length(min_dis),length(n_neighbors),2);
        hops=NaN(length(bin),rep,length(min_dis),length(n_neighbors),3);
    end
    for rp=1:rep
        %if not enough data for umap
        [newtrindx]=check_dataforumap(DATAtrindx{rp},prepDATA.Pos{bnr},...
            minN_umap,IDs,batid,randprop,size(data,1));
        if ~isempty(newtrindx{1})
            DATAtrindx{rp}=newtrindx{1};
        end
        disp(['Rep #' num2str(rp) ' of ' num2str(rep)])
        Train=data(DATAtrindx{rp},:);
        %generate labels
        labs={size(Train,2)};
        for f=1:size(Train,2)
            labs{f}=['f' num2str(f)];
        end
        for md=1:length(min_dis)
            for nn=1:length(n_neighbors)
                disp(['MD: ' num2str(min_dis(md)) ' - NN: ' num2str(n_neighbors(nn))])
                disp(['Running umap...'])
                %run template=train data-set
                [reduction, umap]=run_umap(Train,'parameter_names',labs,...
                    'min_dis',min_dis(md),'n_neighbors',n_neighbors(nn),...
                    'verbose','none');
                fh1=figure('Name','umapPlotCheck','NumberTitle','off');
                plot(reduction(:,1),reduction(:,2),'.','Markersize',5)
                drawnow;
                xlims=xlim;
                ylims=ylim;
                close(figure(fh1));
                %if umap projection is ok
                if max(abs(xlims))<maxproj_umap
                    disp(['Running clustering...'])
                    [idx_d,idx_h]=run_clustering(reduction,db_paras,hdb_paras);
                    clN(bnr,rp,md,nn,1)=length(unique(idx_d));
                    clN(bnr,rp,md,nn,2)=length(unique(idx_h));
                    if plotdata
                        %plot cluster results
                        IDX.d=idx_d;
                        IDX.h=idx_h;
                        plotpos.bnr=bnr;
                        plotpos.md=md;
                        plotpos.nn=nn;
                        plotpos.rep=rp;
                        LIMS.xlims=xlims;
                        LIMS.ylims=ylims;
                        if md==1
                        fhdb=figure('Name',['Control_init: Bin' num2str(bin(plotpos.bnr)) '-dbScan-MD' ...
                            num2str(min_dis(plotpos.md)) '-Rep#' num2str(plotpos.rep)],'NumberTitle','off');
                        fhhdb=figure('Name',['Control_init: Bin' num2str(bin(plotpos.bnr)) '-hdbScan-MD' ...
                            num2str(min_dis(plotpos.md)) '-Rep#' num2str(plotpos.rep)],'NumberTitle','off');
                        end
                        plot_cluster_results(reduction,bios,IDX,DATAtrindx{rp},...
                            min_dis,n_neighbors,plotpos,bin,LIMS,bioflag,fhdb,fhhdb)
                    end
                    disp(['Running hopkins...'])
                    %hopkins cluster evaluation
                    H=NaN(rephop,1);
                    for rh=1:rephop
                        seq=randperm(size(reduction,1));
                        seq=seq(1:round(size(reduction,1)*prop));
                        H(rh)=hopkins_stats(reduction(seq,:),min(reduction(seq,:),[],1), ...
                            max(reduction(seq,:),[],1),round(prop*size(reduction,1))-...
                            round(portion*(prop*size(reduction,1))),distribution,n_neighborhop,flip);
                    end
                    %bootstrap hopkins results and get CI
                    CI=bootci(1e3,@mean,H);
                    hops(bnr,rp,md,nn,:)=[mean(H);CI];
                else
                    disp('umap projection not successful...')
                end
            end
        end
    end
end

