clear all;
close all;

%Example application of how to use the init_umaptemplaterun to identify
%optimal minimal distance and nearest neighbor pair and ideal clustering method;
%and the run_umaptemplate_clustering for comparing data sets
pathname=[pwd '\umap_data\Julie_umap_data_saldeaf_spec_tsc.mat'];
%different min_dis and n_neighbor values
PARAS.min_dis=[0.051 .1 .2];
PARAS.n_neighbors =[3 5 10];
%minimal data for umap
PARAS.minN_umap=100;
%maximal valid projection range for umap
PARAS.maxproj_umap=500;
%to grab all data for calculation set to 0, else 1:
PARAS.randomize_data=1;
%proportion of data to grab randomly if independent of individual
PARAS.randprop=.4;
%randomize weighted among individuals (=same number of calls per
%individual) set to 1, else set 0
PARAS.rand_weighted=1;
%# of repetitions for randomly grabbing N calls for each individual for
%training data-set
PARAS.rep=10;
%distribution of random sample data for hopkins evaluation
PARAS.distribution='uniform_convex_hull';
%neighbor value for hopkins stats
PARAS.n_neighborhop=1;
%flip hopkins
PARAS.flip=0;
%proportion of umap data into hopkins
PARAS.prop=.8;
%proportion of umap data not used for hopkins evaluation
PARAS.portion=0.05;
%# of repetitions for hopkins calculation on percentage (=prop) of umap data
%-prop of that data
PARAS.rephop=20;
%values for clustering
db_paras.md=1;
db_paras.nn=5;
hdb_paras.md=10;
hdb_paras.nn=4;
PARAS.db_paras=db_paras;
PARAS.hdb_paras=hdb_paras;
%plots biosound parameters for each vocalization or snippet (had also to be set
%1 in prep_audio_umap)
bioflag=1;
plotdata=1;
load(pathname)
%run initial umap to find optimal minimal distance and nearest neighbot pair
%and cluster method
[hops,clN]=init_umaptemplaterun(prepDATAs,bin,PARAS,plotdata,bioflag);
save([pathname(1:end-4) '_EVALinit.mat'],'hops','clN','bin')
%plot the hopkin distribution
plot_init_umap_result(pathname,PARAS,bin)
%define optimal minimal distance and nearest neighbot pair from above
%application
PARAS.min_dis=input('Enter optimal minimal distance value: ');
PARAS.n_neighbors=input('Enter optimal nearest neighbor value: ');
PARAS.rep=1;
plotdata=0;
scan='dbscan';% or hdbscan
cmpDATA{1}=prepDATAd;
%run umap with template generation and test data set(s)
[ClusComps,ClusQuals]=run_umaptemplate_clustering(prepDATAs,cmpDATA,bin,PARAS,scan,plotdata,bioflag);
save([pathname(1:end-4) '_EVAL.mat'],'ClusComps','bin')
%plot the cluster distributions for each repetition
%which comparison
compN=1;
plot_cluster_dis(pathname,bin,compN)
%plot comparison of cluster occurrance between saline and deaf for each
%repetition and bin
cluscomp=ClusComps{compN};
clstr={'k';'r';'b';'g';'m';'c';'y';'k'};
fh=figure('Name','Cluster comparison','NumberTitle','off');
for bnr=1:length(bin)
    for rp=1:size(cluscomp,2)        
       figure(fh);
       plot(squeeze(cluscomp(bnr,rp,1,:)),squeeze(cluscomp(bnr,rp,2,:)),...
            [clstr{bnr} '.'],'Markersize',20)
        hold on;
        ylabel('Saline cluster # (%)')
        xlabel('Deaf cluster # (%)')
    end
end
