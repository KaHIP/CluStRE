/******************************************************************************
 * modularity.cpp
 *****************************************************************************/

#include "modularity.h"
#include "partition/onepass_partitioning/floating_block.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"

onepass_modularity::onepass_modularity(PartitionID k0, PartitionID kf,
                                       PartitionID max_blocks, NodeID n_threads,
                                       int mode,
                                       bool hashing /*=false*/, float gamma /*=1.5*/)
    : vertex_partitioning(k0, kf, max_blocks, n_threads, mode, hashing) {
  this->gamma = gamma;
}

onepass_modularity::~onepass_modularity() {}

double onepass_modularity::compute_score(floating_block &block, int my_thread, HeiClus::PartitionConfig &config, PartitionID cluster, int restreaming, LongNodeID curr_node_id, std::vector<floating_block> & blocks) {
  
  double score;

  if(blocks.size() >= config.max_num_clusters && !restreaming) {
    // Compute gain of the modularity when node first initialised to cluster ID : max_num_clusters
    score = block.get_modularity_max_cluster(my_thread, config, current_graph_edge_count, streamed_curr_node_edge_count, curr_node_id, blocks);
  } else {
    if(!restreaming) {
      // Compute the gain in modularity if the streamed node is a singleton and max_num_clusters has not yet been reached.
      score = block.get_modularity_obj(my_thread, config, current_graph_edge_count, streamed_curr_node_edge_count);
    } else {
      // Compute gain of modularity if streamed node is NOT necessarily a singleton (volume of current cluster must be considered)
      // and max_num_clusters has not yet been reached.
      score = block.get_modularity_restream(my_thread, config, current_graph_edge_count, streamed_curr_node_edge_count, curr_node_id, blocks);
    }
  }
  return score;
}

double onepass_modularity::calculate_overall_score(std::vector<NodeID> &blocks_weight, std::vector<std::pair<EdgeWeight,EdgeWeight>> &blocks_sigma, LongEdgeID nmbEdges) {
  LongEdgeID tot_edges = nmbEdges * 2;
  double score = 0.0;
  for(auto &block : blocks_sigma) {
    double first = (static_cast<double>(block.first) / static_cast<double>(tot_edges));
    double second = std::pow((static_cast<double>(block.second) / static_cast<double>(tot_edges)),2);
    score += first - second;
  }

  return score;
}