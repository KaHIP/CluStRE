/******************************************************************************
 * vertex_partitioning.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#include "vertex_partitioning.h"
#include "partition/onepass_partitioning/floating_block.h"

vertex_partitioning::vertex_partitioning(PartitionID k0, PartitionID kf,
                                         PartitionID max_blocks,
                                         NodeID n_threads,
                                         int mode,
                                         bool hashing /*=false*/) {
  this->k0 = k0;
  this->kf = kf;
  this->max_blocks = max_blocks;
  this->n_threads = n_threads;
  this->neighbor_blocks.resize(n_threads);
  original_problem = NULL;
  subproblem_tree_root = NULL;
  parent_block = NULL;
  this->hashing = hashing;
  this->use_self_sorting_array = false;
  int largest_dim = max_blocks;
  EdgeWeight current_graph_edge_count = 0;
  EdgeWeight streamed_curr_node_edge_count = 0;
  this->quotient_edge_count = 0;
  this->overall_delta_mod = 0;
  if(mode != LIGHT_PLUS && mode != LIGHT ) {
    //quotient = new robin_hood::unordered_flat_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash>();
    quotient = new absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash>();

  }
}


double vertex_partitioning::compute_score(floating_block &block, int my_thread, HeiClus::PartitionConfig &config, PartitionID cluster, int restreaming, LongNodeID curr_node_id, std::vector<floating_block> & blocks) {
  return 0;
}

void vertex_partitioning::update_quotient_graph(PartitionID block,
                                                std::vector<std::pair<PartitionID, EdgeWeight>> &neighbor_blocks) {

  // Loop over the neighbours of the streamed node to update / insert edges in the quotient Graph.                              
  for (auto &neighbour_block : neighbor_blocks) {
    if(neighbour_block.first == -1) {
      break;
    }

    PartitionID neighbour = neighbour_block.first;
    EdgeWeight weight = neighbour_block.second;

    std::pair<PartitionID, PartitionID> key = (block <= neighbour) 
             ? std::make_pair(block, neighbour) 
             : std::make_pair(neighbour, block);
             
     // Insert or update the entry in the quotient graph
     EdgeWeight &current_weight = (*quotient)[key];
     if (current_weight == 0 && block != neighbour) {
        // Increment edge count only for new edges (excluding self-loops)
         quotient_edge_count++;
     }
     current_weight += weight;
  }
}
