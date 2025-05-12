/******************************************************************************
 * vertex_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef EXTCLUSTERING_7I4IR31Y
#define EXTCLUSTERING_7I4IR31Y

#include "definitions.h"
#include "random_functions.h"
#include "timer.h"
#include <algorithm>
#include <omp.h>

#include "robin_hood.h"
#include "absl/container/flat_hash_map.h"

#include "data_structure/priority_queues/self_sorting_monotonic_vector.h"
#include "data_structure/hashmap.h"
#include "partition/onepass_partitioning/floating_block.h"

#include "cpi/run_length_compression.hpp"

#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "data_structure/compression_vectors/BatchRunLengthCompression.h"


#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

class extclustering {
public:
    extclustering(LongNodeID rle_length);

    ~extclustering();
    
    virtual void start_clustering(  absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> *& quotient,
                                    std::vector<floating_block> & orig_block,
                                    std::vector <NodeWeight> & stream_blocks_weight,
                                    std::vector <PartitionID> & stream_nodes_assign,
                                    NodeID num_nodes,
                                    EdgeID num_edges,
                                    const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments ) = 0;

    virtual void writeClusteringToDisk(    HeiClus::PartitionConfig &config,
                                           const std::string &filename,
                                           const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments ) = 0;

protected:

    virtual void convert_q_to_local_ds( absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> *& quotient,
                                        std::vector <NodeWeight> & stream_blocks_weight,
                                        NodeID num_nodes,
                                        EdgeID num_edges ) = 0;
    
    virtual void convert_local_ds_to_q( std::vector<floating_block> & orig_block,
                                        std::vector<NodeWeight> & stream_blocks_weight,
                                        std::vector <PartitionID> & stream_nodes_assign,
                                        NodeID num_nodes,
                                        EdgeID num_edges,
                                        const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) = 0;
    
    LongNodeID m_rle_length;
};

#endif /* end of include guard: VERTEX_PARTITIONING_7I4IR31Y */
