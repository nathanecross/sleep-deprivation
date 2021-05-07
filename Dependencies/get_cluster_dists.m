function cluster_dists = get_cluster_dists(Z, num_clusters)
% parse the SAHN tree Z to obtain the distances associated with nodes
% corresponding to clusters upon partitioning the SAHN tree into
% num_clusters clusters: cluster_dists(i) is the distance of the i'th such
% cluster (the numbering of the clusters corresponds to the numbering
% assigned to clusters by the function cluster

num_elements = length(Z) + 1;

cluster_dists = zeros(num_clusters,1);

if(num_clusters == num_elements)
    % all clusters are singletons and hence have 0 associated with their nodes
    % so nothing more to do
   return;
end

pre_sig_node = num_elements - num_clusters;
cur_cluster_num = 1;

for(i=1:(num_clusters-1))
    node1idx = Z(i+pre_sig_node,1);
    node2idx = Z(i+pre_sig_node,2);

    [dist1 clust1idx cur_cluster_num] = get_dist_idcs(node1idx, Z, num_clusters, cur_cluster_num);
    [dist2 clust2idx cur_cluster_num] = get_dist_idcs(node2idx, Z, num_clusters, cur_cluster_num);
    
    cluster_dists(clust1idx) = dist1;
    cluster_dists(clust2idx) = dist2;
end


function [dist_val clust_idx nci] = get_dist_idcs(node_idx, the_SAHN_tree, nc, cci)
% helper function to get distance for cluster and which elements are in
% this cluster

dist_val = -1;
num_elements = size(the_SAHN_tree,1) + 1;
num_nodes_sup = 2*num_elements-(nc-1);

if(node_idx <= num_elements) % cluster is a singleton
    dist_val = 0;
elseif(node_idx < num_nodes_sup)
    dist_val = the_SAHN_tree(node_idx-num_elements,3);
end

if(dist_val < 0) % no cluster at all
    nci = cci;
    clust_idx = [];
else % we have explored a cluster
    clust_idx = cci;
    nci = cci + 1;
end