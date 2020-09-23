function [ClusComps,ClusQuals]=run_umaptemplate_clustering(Template_prepDATA,Test_prepDATA,bin,PARAS,scan,plotdata,bioflag)
%example of supervised umap on time-scaled data prepped through
%prep_audio_umap.m; runs supervised UMAP, template (training data-set) is generated from defined
%number of vocalizations for each saline individual for each existing
%bin and for  fixed min_dis and n_neighbor values retrieved through intial
%clustering testing

min_dis=PARAS.min_dis;
n_neighbors=PARAS.n_neighbors;
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
cluscomp=NaN(length(bin),rep,2,100);
clusqual=NaN(length(bin),rep,2,100,2);
ClusComps={size(Test_prepDATA,2)};
ClusQuals={size(Test_prepDATA,2)};
%run supervised UMAP
for comptest=1:size(Test_prepDATA,2)
    testprepDATA=Test_prepDATA{comptest};
    for bnr=1:length(bin)
        disp(['Bin ' num2str(bin(bnr))])
        datatemp=Template_prepDATA.umap{bnr};
        datatest=testprepDATA.umap{bnr};
        biostest=testprepDATA.bios{bnr};
        biostemp=Template_prepDATA.bios{bnr};
        %generates ID for each snippet
        allindx=NaN(2,5);
        for ids=1:2
            if ids==1
                IDs=Template_prepDATA.IDs;
            else
                IDs=testprepDATA.IDs;
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
                IDs=Template_prepDATA.IDs;
                callN=round(randprop*length(IDs));
            end
            [tempDATAtrindx]=get_random_calls(Template_prepDATA.Pos{bnr},...
                Template_prepDATA.IDs,...
                indxidsall,callN,rep,rand_weighted);
            if rand_weighted~=1
                IDs=testprepDATA.IDs;
                callN=round(randprop*length(IDs));
            end
            [testDATAtrindx]=get_random_calls(testprepDATA.Pos{bnr},...
                testprepDATA.IDs,indxiddall,callN,rep,rand_weighted);
        else
            %grab all data
            rep=1;
            tempDATAtrindx{1}=1:size(datatemp,1);
            testDATAtrindx{1}=1:size(datatest,1);
            cluscomp=NaN(length(bin),rep,2,100);
            clusqual=NaN(length(bin),rep,2,100,2);
        end
        %run the number repetitions
        for rp=1:rep
            
            %if not enough data for umap
            [newtrindx]=check_dataforumap(tempDATAtrindx{rp},Template_prepDATA.Pos{bnr},...
                minN_umap,Template_prepDATA.IDs,indxidsall,randprop,size(datatemp,1));
            if ~isempty(newtrindx{1})
                tempDATAtrindx{rp}=newtrindx{1};
            end
            [newtrindx]=check_dataforumap(testDATAtrindx{rp},testprepDATA.Pos{bnr},...
                minN_umap,testprepDATA.IDs,indxiddall,randprop,size(datatest,1));
            if ~isempty(newtrindx{1})
                testDATAtrindx{rp}=newtrindx{1};
            end
            disp(['Rep #' num2str(rp) ' of ' num2str(rep)])
            Train=datatemp(tempDATAtrindx{rp},:);
            Test=datatest(testDATAtrindx{rp},:);
            %generate labels
            labs={size(Train,2)};
            for f=1:size(Train,2)
                labs{f}=['f' num2str(f)];
            end
            disp(['Running umaps...'])
            %run template=train data-set
            [reduction, umap]=run_umap(Train,'parameter_names',labs,...
                'min_dis',min_dis,'n_neighbors',n_neighbors,...
                'verbose','none');
            fh1=figure('Name','umapPlotCheck','NumberTitle','off');
            plot(reduction(:,1),reduction(:,2),'.','Markersize',5)
            xlims=xlim;
            close(figure(fh1));
            %if umap projection is ok
            if max(abs(xlims))<maxproj_umap
                save([pwd '\umap_data\UMAP_saltrain_tsc.mat'], 'umap','reduction','-v7.3');
                %run test with template=train data-set
                [reductiont, umapt]= run_umap(Test, 'template_file',...
                    [pwd '\umap_data\UMAP_saltrain_tsc.mat'],'parameter_names',labs,...
                    'verbose','none');
                disp(['Running clustering...'])
                vec=[reduction;reductiont];
                clindx{1}=1:size(reduction,1);
                clindx{2}=size(reduction,1)+1:size(vec,1);
                if strcmp(scan,'dbscan')
                    %dbscan
                    D=pdist2(vec,vec);
                    [idx, corepts] = dbscan(D,1,5,'Distance','precomputed');
                elseif strcmp(scan,'hdbscan')
                    %hdbscan
                    clusteringh=HDBSCAN(vec);
                    clusteringh.minpts=10;
                    clusteringh.minclustsize=4;
                    fit_model(clusteringh)
                    get_best_clusters(clusteringh)
                    get_membership(clusteringh)
                    idx=double(clusteringh.labels);
                end
                for mi=1:2
                    if max(idx)>0
                        [nbins,idxh]=hist(idx(clindx{mi}),-1:1:max(idx));
                        cluscomp(bnr,rp,mi,1:length(nbins))=nbins./length(clindx{mi});
                        minmxidx=[min(idx) max(idx)];
                        %assess cluster qualities
                        for cln=minmxidx(1): minmxidx(2)
                            idics=idx(clindx{mi});
                            clnindx=find(idics==cln);
                            if length(clnindx)>2
                                [CluSep, ~] = Cluster_Quality(vec,idics(clnindx));
                                clusqual(bnr,rp,mi,cln,:)=[CluSep.IsolationDist CluSep.L];
                            end
                        end
                        if plotdata
                            %plot umap reduction, template and test
                            if mi==1
                                fh1=figure('Name','umapPlotCheck','NumberTitle','off');
                                plot(vec(:,1),vec(:,2),'.','Markersize',5)
                                title(['Bin :' num2str(bin(bnr))])
                                xlims=xlim;
                                ylims=ylim;
                                close(figure(fh1))
                            end
                            %plot clustering of umap reduction
                            plotpos.mi=mi;
                            plotpos.bnr=bnr;
                            plotpos.compN=comptest;
                            plotpos.rep=rp;
                            LIMS.xlims=xlims;
                            LIMS.ylims=ylims;
                            if mi==1
                                fhclus=figure('Name',['Cluster: Bin' num2str(bin(bnr)) '-Rep#' num2str(plotpos.rep)],'NumberTitle','off');
                                plot_umap_template_test(vec,Train,idx,clindx{mi},bin,...
                                    minmxidx,0,[biostemp;biostest],...
                                    biostemp(tempDATAtrindx{rp},:),Template_prepDATA.fhigh,...
                                    plotpos,LIMS,bioflag,fhclus)
                            else
                                plot_umap_template_test(vec,Test,idx,clindx{mi},bin,...
                                    minmxidx,length(clindx{1}),[biostemp;biostest],...
                                    biostest(testDATAtrindx{rp},:),testprepDATA.fhigh,...
                                    plotpos,LIMS,bioflag,fhclus)
                            end
                        end
                    else
                        disp('Clustering not successful...')
                    end
                end
            else
                disp('umap projection not successful...')
            end
        end
    end
    ClusComps{comptest}= cluscomp;
    ClusQuals{comptest}=clusqual;
end