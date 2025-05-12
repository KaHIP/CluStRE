/******************************************************************************
 * modularity.h
 *****************************************************************************/

#ifndef ONEPASS_MODULARITY_7I4IR31Y
#define ONEPASS_MODULARITY_7I4IR31Y

#include <algorithm>

#include "definitions.h"
#include "random_functions.h"
#include "timer.h"

#include "vertex_partitioning.h"

class onepass_modularity : public vertex_partitioning {
public:
    onepass_modularity(PartitionID k0, PartitionID kf, PartitionID max_blocks,
                       NodeID n_threads, int mode, bool hashing = false, float gamma = 1.5);

    virtual ~onepass_modularity();

    double calculate_overall_score(std::vector<NodeID> &blocks_weight, std::vector<std::pair<EdgeWeight,EdgeWeight>> &blocks_sigma, LongEdgeID nmbEdges);

protected:
    double compute_score(floating_block &block, int my_thread, HeiClus::PartitionConfig &config, PartitionID cluster, int restreaming, LongNodeID curr_node_id, std::vector<floating_block> & blocks);
    float gamma;
};

#endif /* end of include guard: ONEPASS_MODULARITY_7I4IR31Y */
