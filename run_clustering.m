function [idx_d,idx_h]=run_clustering(DATA,db_paras,hdb_paras)

%dbscan
D=pdist2(DATA,DATA);
[idx_d, corepts] = dbscan(D,db_paras.md,db_paras.nn,'Distance','precomputed');
%hdbscan
clusteringh=HDBSCAN(DATA);
clusteringh.minpts=hdb_paras.md;
clusteringh.minclustsize=hdb_paras.nn;
fit_model(clusteringh)
get_best_clusters(clusteringh)
get_membership(clusteringh)
idx_h=clusteringh.labels;
