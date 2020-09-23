clear all;
close all;
%example of supervised umap on time-scaled data prepped through
%prep_audio_umap.m; runs supervised UMAP, template (training data-set)is generated from defined
%number of vocalizations for specific saline and deaf individuals for each existing
%bin and same for the test data-set

addpath([pwd '\umap\'])
addpath([pwd '\umap\sgdCpp_files'])
addpath([pwd '\util\'])
%generate various combinations of train bat ID numbers and corresponding
%test ID numbers for both saline and deaf (5 bats each)
trids=[1 2 3;2 3 4;3 4 5;1 2 4;1 2 5;1 3 4;1 3 5;1 4 5;2 4 5;2 3 5];
tsids=[4 5;  1 5;  1 2;   3 5; 3 4;  2 5;  2 4;   2 3;  1 3;  1 4];
%# of repetitions for randomly grabbing N calls for each individual for both
%training and test data-set
rep=50;
%plot results (=1)
plotdata=0;
%calculate the probablity density function between umap projections (=1)
%for now this is instead of using labels outputted by umap
runpdf=1;
%example of prepped data
load([pwd '\umap_data\umap_data_saldeaf_spec_tsc.mat'])
%run supervised UMAP
for bnr=1:length(bin)
    datadeaf=prepDATAd.umap{bnr};
    datasal=prepDATAs.umap{bnr};
    biossal=prepDATAs.bios{bnr};
    biosdeaf=prepDATAd.bios{bnr};
    %generates ID for each snippet
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
        allindx(ids,:)=indxid;
    end
    %define maximal possible number of calls for each individual
    callN=min(min(allindx));
    %run through the different train/test combinations
    for ivar=1:size(trids,1)
        %grab the random calls
        [sDATAtrindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,indxidsall(trids(ivar,:)),callN,rep);
        [sDATAtsindx]=get_random_calls(prepDATAs.Pos{bnr},prepDATAs.IDs,indxidsall(tsids(ivar,:)),callN,rep);
        [dDATAtrindx]=get_random_calls(prepDATAd.Pos{bnr},prepDATAd.IDs,indxiddall(trids(ivar,:)),callN,rep);
        [dDATAtsindx]=get_random_calls(prepDATAd.Pos{bnr},prepDATAd.IDs,indxiddall(tsids(ivar,:)),callN,rep);
        %run the number repetitions for a set
        for rp=1:rep
            Train=[datasal(sDATAtrindx{rp},:);datadeaf(dDATAtrindx{rp},:)];
            Train(1:size(datasal(sDATAtrindx{rp},:),1),end+1)=0;
            Train(size(datasal(sDATAtrindx{rp},:),1)+1:end,end)=1;
            indxtr{1}=1:size(datasal(sDATAtrindx{rp},:),1);
            indxtr{2}=size(datasal(sDATAtrindx{rp},:),1)+1:size(Train,1);
            Test=[datasal(sDATAtsindx{rp},:);datadeaf(dDATAtsindx{rp},:)];
            indxts{1}=1:size(datasal(sDATAtsindx{rp},:),1);
            indxts{2}=size(datasal(sDATAtsindx{rp},:),1)+1:size(Test,1);
            %generate labels
            if rp==1 && ivar==1
                labs=[];
                for f=1:size(Train,2)-1
                    labs{f}=['f' num2str(f)];
                end
            end
            %run template=train data-set
            [reduction, umap,clustering]=run_umap(Train,'parameter_names',labs,...
                'label_column',size(Train,2),...
                'qf_dissimilarity',true,'verbose','none');
            save([pwd '\umap_data\UMAP_saldeaftrain_tsc.mat'], 'umap','reduction','-v7.3');
            if plotdata
                %plot the projection for saline and deaf train data-sets
                figure;hold on;colorstr={'k';'r'};tits={'Saline-Deaf-Train'};
                for pl=1:2
                    plot(reduction(indxtr{pl},1),reduction(indxtr{pl},2),[colorstr{pl} '.'],'Markersize',5)
                end
                xlims=xlim;
                ylims=ylim;
                set(gca,'ylim',[ylims],'xlim',[xlims])
                title(['Bin' num2str(bin(bnr)) ', rep' num2str(rp) ': ' tits])
            end
            %run test with template=train data-set
            [reductiont, umapt,clustering]= run_umap(Test, 'template_file',...
                [pwd '\umap_data\UMAP_saldeaftrain_tsc.mat'],'parameter_names',labs,...
                'qf_dissimilarity',true,'verbose','none');
            %  save(['C:\Users\Daria\audio_stuff\umap_data\UMAP_saldeaftest_' num2str(bin(bnr)) '_rep' num2str(rp) '_tsc.mat'],...
            %'umapt','reductiont','xlims','ylims','-v7.3');
            %for now calculate probablity density functions between train
            %sets and test sets to evaluate how well 'classification' of
            %data possible
            if runpdf
                N=20;
                [rx,ry,rpdfstr]=jointpdf(reduction(indxtr{1},1),reduction(indxtr{1},2),N,plotdata);
                [rx,ry,rpdfdtr]=jointpdf(reduction(indxtr{2},1),reduction(indxtr{2},2),N,plotdata);
                [rx,ry,rpdfsts]=jointpdf(reductiont(indxts{1},1),reductiont(indxts{1},2),N,plotdata);
                [rx,ry,rpdfdts]=jointpdf(reductiont(indxts{2},1),reductiont(indxts{2},2),N,plotdata);
                pdfdifftr=(rpdfstr)-(rpdfdtr);
                pdfclass=NaN(size(pdfdifftr));
                pdfclass(find(pdfdifftr<0))=2;
                pdfclass(find(pdfdifftr>0))=1;
                pdfdiffts=(rpdfsts)-(rpdfdts);
                pdfclasstest=NaN(size(pdfdiffts));
                pdfclasstest(find(pdfdiffts<0))=2;
                pdfclasstest(find(pdfdiffts>0))=1;
                dvals(bnr,ivar,rp)=length(find((pdfclass-pdfclasstest)==0))/length(find(~isnan(pdfclass-pdfclasstest)))*100;
                if plotdata
                    figure;surf(rx,ry,pdfdiffs)
                    figure;surf(rx,ry,pdfdiffs)
                    figure;hold on;colorstr={'k';'r'};tits={'Saline-Deaf-Test'};
                    for pl=1:2
                        plot(reductiont(indxts{pl},1),reductiont(indxts{pl},2),[colorstr{pl} '.'],'Markersize',5)
                        set(gca,'ylim',[ylims],'xlim',[xlims])
                        title(['Bin' num2str(bin(bnr)) ', rep' num2str(rp) ': ' tits])
                    end
                end
                close all;
            end
        end
    end
end
%save('C:\Users\Daria\umap_stuff\umap_supervised_perc.mat','dvals')
if plotdata
    figure;
    for bnr=1:size(dvals,1)
        pcs=squeeze(dvals(bnr,:,:));
        [nbins,indx]=hist(pcs(:),0:10:100);
        subplot(3,2,bnr);bar(indx,nbins./sum(nbins))
    end
end