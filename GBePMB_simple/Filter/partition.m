function [partitions,clusters,np,n_clusters,partitions_indices] = partition(W,mbm,paras)
%

n_track = length(mbm.track);
%use DBSCAN to obtain multiple partitions
partitions = gen_partitions(W,paras.dbscan,n_track);
%number of partitions
np = length(partitions);

%find all the unique clusters in all measurement partitions
[clusters,IA,IC] = unique(cell2mat(partitions')','rows');
clusters = clusters';
n_clusters = length(IA);
%reconstruct partitions to let it contain indices of clusters
nc_p = cellfun(@(x) size(x,2),partitions);
partitions_indices = cell(np,1);
idx = 0;
for i = 1:np
    partitions_indices{i} = IC(idx+1:idx+nc_p(i));
    idx = idx + nc_p(i);
end

end

