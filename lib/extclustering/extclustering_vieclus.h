/******************************************************************************
 * vertex_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef EXTCLUSTERING_VIECLUS_7I4IR31Y
#define EXTCLUSTERING_VIECLUS_7I4IR31Y

#include "definitions.h"
#include "random_functions.h"
#include "timer.h"
#include <algorithm>
#include <omp.h>

#include "data_structure/priority_queues/self_sorting_monotonic_vector.h"
#include "data_structure/hashmap.h"
#include "partition/onepass_partitioning/floating_block.h"
#include "extclustering/extclustering.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/VieClus/lib/parallel_mh_clustering/parallel_mh_async_clustering.h"
#include "extern/KaHIP/lib/io/graph_io.h"
#include <mpi.h>
#include "absl/container/flat_hash_map.h"


#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

class extclustering_vieclus : public extclustering {
public:
    extclustering_vieclus(LongNodeID rle_length, int argn, char **argv, HeiClus::PartitionConfig & HeiClus_Partconfig);

    ~extclustering_vieclus();
    
    virtual void start_clustering(  absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> *& quotient,
                                    std::vector<floating_block> & orig_block,
                                    std::vector <NodeWeight> & stream_blocks_weight,
                                    std::vector <PartitionID> & stream_nodes_assign,
                                    NodeID num_nodes,
                                    EdgeID num_edges,
                                    const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments );
    
    virtual void writeClusteringToDisk(HeiClus::PartitionConfig &config,
                                           const std::string &filename,
                                           const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);


protected:

    virtual void convert_q_to_local_ds( absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> *& quotient,
                                        std::vector <NodeWeight> & stream_blocks_weight,
                                        NodeID num_nodes,
                                        EdgeID num_edges );
    
    virtual void convert_local_ds_to_q( std::vector<floating_block> & orig_block,
                                        std::vector<NodeWeight> & stream_blocks_weight,
                                        std::vector <PartitionID> & stream_nodes_assign,
                                        NodeID num_nodes,
                                        EdgeID num_edges,
                                        const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);
    
    KaHIP::graph_access * G;
    KaHIP::PartitionConfig * KaHIP_partition_config;
    int argn;
    char **argv;

};

#endif /* end of include guard: VERTEX_PARTITIONING_7I4IR31Y */
