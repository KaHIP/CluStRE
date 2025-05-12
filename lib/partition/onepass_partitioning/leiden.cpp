/******************************************************************************
 * leiden.cpp
 *****************************************************************************/

#include "leiden.h"
#include "partition/onepass_partitioning/floating_block.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"

onepass_leiden::onepass_leiden(PartitionID k0, PartitionID kf,
                               PartitionID max_blocks, NodeID n_threads,
                               int mode,
                               bool hashing /*=false*/, float gamma /*=1.5*/)
    : vertex_partitioning(k0, kf, max_blocks, n_threads, mode, hashing) {
  this->cpm_gamma = gamma;
}

onepass_leiden::~onepass_leiden() {}

double onepass_leiden::compute_score(floating_block &block, int my_thread, HeiClus::PartitionConfig &config, PartitionID cluster, int restreaming, LongNodeID curr_node_id, std::vector<floating_block> & blocks) {
  return block.get_leiden_obj(my_thread, cpm_gamma, config, cluster);
}

double onepass_leiden::calculate_overall_score(std::vector<NodeID> &blocks_weight, std::vector<std::pair<EdgeWeight,EdgeWeight>> &blocks_sigma, EdgeID nmbEdges) {
  float score = 0;

  for(PartitionID i = 0; i < blocks_sigma.size(); i++) {
    score += (static_cast<double>(blocks_sigma[i].first) / 2.0) - cpm_gamma * (static_cast<double>(blocks_weight[i] * (blocks_weight[i] - 1)) / 2.0);
  }

  score *= 2;

  return score;
}

