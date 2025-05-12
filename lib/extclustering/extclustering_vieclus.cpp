/******************************************************************************
 * vertex_partitioning.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/
#include "extclustering/extclustering_vieclus.h"
#include "extclustering/extclustering.h"
#include <queue>

extclustering_vieclus::extclustering_vieclus(LongNodeID rle_length, int argn, char **argv, HeiClus::PartitionConfig & HeiClus_Partconfig)
    : extclustering(rle_length) {

    // Initialise the partition configuration for VieClus 
    
    this->argn = argn;
    this->argv = argv;

    KaHIP_partition_config = new KaHIP::PartitionConfig();
    G = NULL;

    VieClus::configuration cfg;
    cfg.standard((*KaHIP_partition_config));
    cfg.strong((*KaHIP_partition_config));

    KaHIP_partition_config->time_limit = HeiClus_Partconfig.ext_algorithm_time;
    KaHIP_partition_config->cluster_upperbound = std::numeric_limits< NodeWeight >::max()/2;
    KaHIP_partition_config->max_cluster_edge_volume = HeiClus_Partconfig.max_cluster_edge_volume;
}

extclustering_vieclus::~extclustering_vieclus() {}

void extclustering_vieclus::start_clustering(   absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> *& quotient,
                                                std::vector<floating_block> & orig_block,
                                                std::vector<NodeWeight> & stream_blocks_weight,
                                                std::vector <PartitionID> & stream_nodes_assign,
                                                NodeID num_nodes,
                                                EdgeID num_edges,
                                                const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) {
    MPI_Init(&argn, &argv);

    convert_q_to_local_ds(quotient, stream_blocks_weight, num_nodes, num_edges);

    G->set_partition_count(G->get_partition_count_compute());

    omp_set_num_threads(1);

    KaHIP_partition_config->k = 1;
    parallel_mh_async_clustering mh;
    mh.perform_partitioning((*KaHIP_partition_config), (*G));

    int rank, size;
    MPI_Comm communicator = MPI_COMM_WORLD;
    MPI_Comm_rank( communicator, &rank);
    MPI_Comm_size( communicator, &size);

    G->set_partition_count(G->get_partition_count_compute());

    convert_local_ds_to_q(orig_block, stream_blocks_weight, stream_nodes_assign, num_nodes, num_edges, block_assignments);

    MPI_Finalize();

    delete G;
    delete KaHIP_partition_config;
}

void extclustering_vieclus::convert_q_to_local_ds(  absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> *& quotient,
                                                    std::vector<NodeWeight> &stream_blocks_weight,
                                                    NodeID num_nodes,
                                                    EdgeID num_edges) {

    long nmbNodes = num_nodes;
    long nmbEdges = num_edges;

    if (2 * num_edges > std::numeric_limits<int>::max() || num_nodes > std::numeric_limits<int>::max()) {
        std::cerr << "The graph is too large. Currently only 32bit supported!" << std::endl;
        exit(0);
    }

    nmbEdges *= 2; // since we have forward and backward edges

    std::vector<std::vector<std::pair<NodeID, EdgeWeight>>> Q(nmbNodes);

    // Move the upper half of the quotient graph to the static data structure Q,
    // after each inserted edge, delete it from the quotient.
    // Q is built, at the same time the dynamic quotient data structure is deleted.
    for (auto it = quotient->begin(); it != quotient->end(); ) {
        Q[it->first.first].push_back({it->first.second, it->second});
        //it = quotient->erase(it);
        auto current = it++;
        quotient->erase(current);
    }

    delete quotient;

    // Create Static Data structrue of VieClus
    G = new KaHIP::graph_access();

    G->start_construction(nmbNodes, nmbEdges);
    G->resizeSelfLoops(nmbNodes);

    int counter = 0;

    // Initialise the VieClus data structure. At the same time, Clear the static data structure Q.
    for(PartitionID s_cluster = 0; s_cluster < Q.size(); s_cluster++) {
        NodeID node = G->new_node();
        G->setPartitionIndex(node, counter++);
        G->setNodeWeight(node, stream_blocks_weight[node]);

        for (PartitionID t_cluster = 0; t_cluster < Q[s_cluster].size(); t_cluster++) {
            NodeID target = Q[s_cluster][t_cluster].first;
            EdgeWeight edge_weight = Q[s_cluster][t_cluster].second;

            // VieClus needs both edges ( 1 -> 2 AND 2 -> 1 )
            if(node < target) {
                Q[target].push_back({node, edge_weight});
            }

            if(target == node) {
                G->setSelfLoop(target, edge_weight * 2);
                continue;
            }
            EdgeID e = G->new_edge(node, target);
            G->setEdgeWeight(e, edge_weight);
        }
        
        //  Swap the contents of Q[s_cluster] with an empty temporary (tmp). 
        //  After the swap, Q[s_cluster] will be empty with zero capacity, and tmp will contain the old data 
        //  (which will then be out of scope and should free its allocation).
        std::vector<std::pair<NodeID, EdgeWeight>>().swap(Q[s_cluster]);
    }
    
    G->finish_construction();
}

void extclustering_vieclus::convert_local_ds_to_q(  std::vector<floating_block> & orig_block,
                                                    std::vector<NodeWeight> & stream_blocks_weight,
                                                    std::vector <PartitionID> & stream_nodes_assign,
                                                    NodeID num_nodes,
                                                    EdgeID num_edges,
                                                    const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) {
    
    if(m_rle_length == -1) {
        for(NodeID node = 0; node < stream_nodes_assign.size(); node++) {
            stream_nodes_assign[node] = G->getPartitionIndex(stream_nodes_assign[node]);
        }
    }

    // Update block weights vector to correspond to the amount of clusteres outputed by VieClus
    while(stream_blocks_weight.size() > G->get_partition_count()) {
        stream_blocks_weight.pop_back();
    }

    for(int i = 0; i < stream_blocks_weight.size(); i++) {
        stream_blocks_weight[i] = 0;
    }

    // Update Volume of Clusters
    for(NodeID node = 0; node < G->number_of_nodes(); node++) {
        stream_blocks_weight[G->getPartitionIndex(node)] += G->getNodeWeight(node);
        if(node != G->getPartitionIndex(node)) {
            orig_block[G->getPartitionIndex(node)].sigma_cluster += orig_block[node].sigma_cluster;
            orig_block[node].sigma_cluster = 0;
        }
    }
}

void extclustering_vieclus::writeClusteringToDisk(HeiClus::PartitionConfig &config,
                                           const std::string &filename,
                                           const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) {
    
    std::ofstream f(filename.c_str());
    //std::cout << "writing cluster to " << filename << " ... " << std::endl;

    if(config.rle_length == 0) {
        for (int node = 0; node < config.total_nodes; node++) {
            f << G->getPartitionIndex(block_assignments->GetValueByIndex(node)) << "\n";
        }
    } 
}

