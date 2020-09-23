clear all;
close all;
%example of supervised umap on time-scaled data prepped through
%prep_audio_umap.m; runs supervised UMAP, template (training data-set)is generated from defined
%number of vocalizations for each saline individual for each existing
%bin and for different min_dis and n_neighbor values

addpath([pwd '\umap\'])
addpath([pwd '\umap\sgdCpp_files'])
addpath([pwd '\util\'])
%different min_dis and n_neighbor values
min_dis=[0.051 .1 .2 .4 .6 .79];
n_neighbors =[3 5 10 20 50 100 150 199];
%# of repetitions for randomly grabbing N calls for each individual for
%training data-set
rep=1;
%plot results (=1)
plotdata=0;
%example of prepped data
load([pwd '\umap_data\umap_data_saldeaf_spec_tsc.mat'])
%run supervised UMAP
for bnr=1:length(bin)
    datadeaf=prepDATAd.umap{bnr};
    datasal=prepDATAs.umap{bnr};
    biossal=prepDATAs.bios{bnr};
    %generates ID for each snippet
    for ids=1
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
        allindx(ids,:)=indxid;
    end
    %define maximal possible number of calls for each individual
    callN=min(min(allindx));
    %run through the different train/test combinations
    for ivar=1
        %grab the random calls
        [sDATAtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,indxidsall,callN,rep);
        %run the number repetitions for a set
        for rp=1:rep
            Train=datasal(sDATAtrindx{rp},:);
            %generate labels
            if rp==1 && ivar==1
                labs=[];
                for f=1:size(Train,2)
                    labs{f}=['f' num2str(f)];
                end
            end
            for md=1:length(min_dis)
                for nn=1:length(n_neighbors)
                    %run template=train data-set
                    [reduction, umap,clustering]=run_umap(Train,'parameter_names',labs,...
                        'min_dis',min_dis(md),'n_neighbors',n_neighbors(nn),...
                        'verbose','none');
                    clN(md,nn)=length(unique(clustering));
%                     figure(md);subplot(2,4,nn);plot(reduction(:,1),reduction(:,2),...
%                         '.','color','k','Markersize',5)
%                     title(['MD:' num2str(min_dis(md)) '_NN:' num2str(n_neighbors(nn))])
                end
            end
            
            % save([pwd '\umap_data\UMAP_saldeaftrain_tsc.mat'], 'umap','reduction','-v7.3');
        end
    end
end
