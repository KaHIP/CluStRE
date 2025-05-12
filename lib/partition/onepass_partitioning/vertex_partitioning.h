/******************************************************************************
 * vertex_partitioning.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef VERTEX_PARTITIONING_7I4IR31Y
#define VERTEX_PARTITIONING_7I4IR31Y

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

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

class vertex_partitioning {
public:
  vertex_partitioning(PartitionID k0, PartitionID kf, PartitionID max_blocks,
                      NodeID n_threads, int mode, bool hashing = false);
  ~vertex_partitioning() {}
  void instantiate_blocks(LongNodeID n, LongEdgeID m, PartitionID k,
                                  ImbalanceType epsilon);
  void load_edge(PartitionID block, EdgeWeight e_weight, int my_thread);
  void clear_edgeweight(PartitionID block, int my_thread);
  void propagate_load_edge(EdgeWeight e_weight, PartitionID block,
                           int my_thread);
  void propagate_clear_edgeweight(int my_thread);
  PartitionID solve_node(LongNodeID curr_node_id, NodeWeight curr_node_weight,
                         PartitionID previous_assignment, double kappa,
                         int my_thread);
  PartitionID set_decision(PartitionID block, LongNodeID curr_node_id,
                           NodeWeight curr_node_weight, int my_thread, LongEdgeID degree);
  void set_original_problem(vertex_partitioning *original_problem);
  void set_parent_block(floating_block *parent_block);
  void clear_edgeweight_blocks(
      std::vector<std::pair<PartitionID, EdgeWeight>> &neighbor_blocks_thread,
      PartitionID n_elements, int my_thread);

  // methods for dealing with blocks
  void load_edge_block(floating_block &real_block, EdgeWeight e_weight,
                       int my_thread);
  void clear_edgeweight_block(floating_block &real_block, int my_thread);
  void assign_node_block(floating_block &real_block,
                                LongNodeID curr_node_id,
                                NodeWeight curr_node_weight, int my_thread, LongEdgeID degree);
  PartitionID force_decision_block(floating_block &real_block,
                                   LongNodeID curr_node_id,
                                   NodeWeight curr_node_weight, int my_thread);
  floating_block *get_block_address(PartitionID block_id) const;
  void enable_self_sorting_array();

  //Streaming Graph Clustering
  void update_quotient_graph(PartitionID block, 
                              std::vector<std::pair<PartitionID, EdgeWeight>> &neighbor_blocks);
  
  PartitionID solve_node_clustering(LongNodeID curr_node_id, NodeWeight curr_node_weight, int restreaming,
                         int my_thread, std::vector<std::pair<PartitionID, EdgeWeight>> &neighbours,
                         std::vector<PartitionID> &clusters_to_ix_mapping,
                         HeiClus::PartitionConfig &config,
                         PartitionID previous_assignment, double kappa, bool local_search);
  
  void increment_graph_edge_count(EdgeWeight current_graph_edge_count, int restreaming);
  void increment_curr_node_edge(EdgeWeight node_edge_count);
  void increment_cluster_sigma(PartitionID clusterID);
  void reset_streamed_edge_count();

  virtual double calculate_overall_score(std::vector<NodeID> &blocks_weight, std::vector<std::pair<EdgeWeight,EdgeWeight>> &blocks_sigma, LongEdgeID nmbEdges);

  double overall_delta_mod;
  double curr_round_delta_mod;
  std::vector<floating_block> blocks;
  //robin_hood::unordered_flat_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> * quotient = NULL;
  absl::flat_hash_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> * quotient = NULL;
  EdgeID quotient_edge_count;

protected:

  PartitionID solve_linear_complexity(LongNodeID curr_node_id,
                                      NodeWeight curr_node_weight,
                                      int my_thread,
                                      HeiClus::PartitionConfig &config);
  virtual double compute_score(floating_block &block, int my_thread, HeiClus::PartitionConfig &config, PartitionID cluster, int restreaming, LongNodeID curr_node_id, std::vector<floating_block> & blocks);

  //Streaming Graph Clustering
  PartitionID solveClustering(LongNodeID curr_node_id,
                              NodeWeight curr_node_weight, 
                              int restreaming,
                              int my_thread,
                              std::vector<std::pair<PartitionID,EdgeWeight>> &neighbours,
                              std::vector<PartitionID> &clusters_to_ix_mapping,
                              HeiClus::PartitionConfig &config,
                              PartitionID previous_assignment, double kappa, bool local_search); 

  NodeID n_threads;
  PartitionID k0;
  PartitionID kf;
  PartitionID max_blocks;
  ImbalanceType base_size_constraint;
  bool hashing;
  std::vector<std::vector<PartitionID>> neighbor_blocks;
  HeiClus::random_functions::fastRandBool<uint64_t> random_obj;
  self_sorting_monotonic_vector<PartitionID, NodeWeight> sorted_blocks;
  bool use_self_sorting_array;

  vertex_partitioning *original_problem;
  vertex_partitioning *subproblem_tree_root;
  floating_block *parent_block;

  //streaming graph clustering
  EdgeWeight current_graph_edge_count;
  EdgeWeight streamed_curr_node_edge_count;
};

inline void vertex_partitioning::clear_edgeweight_blocks(
    std::vector<std::pair<PartitionID, EdgeWeight>> &neighbor_blocks_thread,
    PartitionID n_elements, int my_thread) {
  for (PartitionID key = 0; key < n_elements; key++) {
    auto &element = neighbor_blocks_thread[key];
    clear_edgeweight(element.first, my_thread);
  }
}

inline void vertex_partitioning::load_edge(PartitionID block,
                                           EdgeWeight e_weight, int my_thread) {
  load_edge_block(blocks[block], e_weight, my_thread);
}

inline void vertex_partitioning::load_edge_block(floating_block &real_block,
                                                 EdgeWeight e_weight,
                                                 int my_thread) {
  real_block.increment_edgeweight(my_thread, e_weight);
  //real_block.parent_problem->propagate_load_edge(
  //    e_weight, real_block.get_block_id(), my_thread);
}

inline void vertex_partitioning::propagate_load_edge(EdgeWeight e_weight,
                                                     PartitionID block,
                                                     int my_thread) {
  neighbor_blocks[my_thread].push_back(block);
}

inline void vertex_partitioning::clear_edgeweight(PartitionID block,
                                                  int my_thread) {
  clear_edgeweight_block(blocks[block], my_thread);
}

inline void
vertex_partitioning::clear_edgeweight_block(floating_block &real_block,
                                            int my_thread) {
  real_block.set_edgeweight(my_thread, 0);
  //real_block.parent_problem->propagate_clear_edgeweight(my_thread);
}

inline void vertex_partitioning::propagate_clear_edgeweight(int my_thread) {
  neighbor_blocks[my_thread].clear();
}

inline PartitionID
vertex_partitioning::set_decision(PartitionID block, LongNodeID curr_node_id,
                                  NodeWeight curr_node_weight, int my_thread, LongEdgeID degree) {
  assign_node_block(blocks[block], curr_node_id, curr_node_weight,
                           my_thread, degree);
  return;
}

inline void vertex_partitioning::assign_node_block(
    floating_block & real_block, LongNodeID curr_node_id,
    NodeWeight curr_node_weight, int my_thread, LongEdgeID degree) {

  //PartitionID decision = real_block.get_block_id();
  //real_block.increment_curr_weight(curr_node_weight);
  real_block.update_cluster_sigma(degree);
  //return decision;
}

inline PartitionID vertex_partitioning::solve_node_clustering(LongNodeID curr_node_id, 
                                                              NodeWeight curr_node_weight,
                                                              int restreaming,
                                                              int my_thread, 
                                                              std::vector<std::pair<PartitionID, EdgeWeight>> &neighbours,
                                                              std::vector<PartitionID> &clusters_to_ix_mapping,
                                                              HeiClus::PartitionConfig &config,
                                                              PartitionID previous_assignment, 
                                                              double kappa,
                                                              bool local_search) {
  
  PartitionID decision;
  decision = solveClustering(curr_node_id, curr_node_weight, restreaming, my_thread, neighbours, clusters_to_ix_mapping, config, previous_assignment, kappa, local_search);

  return decision;
}

inline PartitionID vertex_partitioning::solveClustering(LongNodeID curr_node_id,
                                                        NodeWeight curr_node_weight,
                                                        int restreaming,
                                                        int my_thread,
                                                        std::vector<std::pair<PartitionID,EdgeWeight>> &neighbours,
                                                        std::vector<PartitionID> &clusters_to_ix_mapping,
                                                        HeiClus::PartitionConfig &config,
                                                        PartitionID previous_assignment, 
                                                        double kappa,
                                                        bool local_search) {
  double best;
  double curr_score;
  PartitionID decision;

  // Set the scores and initial block assignments to streamed vector aligned with current graph setting
  if(blocks.size() >= config.max_num_clusters && !restreaming) {
    best = -2;
    curr_score = -2;
    decision = blocks.size()-1;
  } else if (restreaming && blocks.size() >= config.max_num_clusters) {
    best = 0;
    curr_score = 0;
    decision = (*config.stream_nodes_assign)[curr_node_id];
  } else if (restreaming) {
    best = 0;
    curr_score = 0;
    decision = (*config.stream_nodes_assign)[curr_node_id];
  } else {
    best = 0;
    curr_score = 0;
    decision = blocks.size();
  }
 
  // Loop over neighbours of streamed node
  for (auto &block : neighbours) {
	  if(block.first == -1) {
      		break;
    }

    // Respect edge Volume limit if specified
    if (blocks[block.first].sigma_cluster + streamed_curr_node_edge_count > config.max_cluster_edge_volume) {
      continue;
    }
    
    // streamed node and neighbour (block.first) have the same cluster id. no need to calculate gain.
    if(restreaming && (*config.stream_nodes_assign)[curr_node_id] == block.first) {
      continue;
    }
    
    // Compute modularity gain when assigning a node the same cluster ID as its neighbour (block.first)
    curr_score = compute_score(blocks[block.first], my_thread, config, block.first, restreaming, curr_node_id, blocks);
    if (block.first == previous_assignment) {
      // Boost previously assigned cluster if kappa is initialised kappa != 1.
      if (curr_score > 0) {
        curr_score = curr_score * kappa;
      } else {
        curr_score = curr_score / kappa;
      }
    }
    if (curr_score > best) {
      decision = block.first;
      best = curr_score;
    }
  }

  // If the streamed node has no neighbours and the maximum number of clusters has been reached, 
  // assign the streamed vertex to the cluster ID: max_num_clusters.
  if(best <= -2 && blocks.size() >= config.max_num_clusters) {
	  best = 0;
	  decision = blocks.size()-1;
  }
  
  // Update the overall estimated modularity score, which will then be used for the stopping criteria.
  // in the local search phase.
  if(!local_search) {
    overall_delta_mod += best;
  } else {
    curr_round_delta_mod += best;
  }

  // Streamed node has been assigned to its own cluster
  if (decision == blocks.size()) {
    blocks.emplace_back(std::move(floating_block(this, blocks.size(), n_threads)));

    std::pair<PartitionID, EdgeWeight> new_neighbour_block = {-1, 0};
    config.neighbor_blocks[my_thread].emplace_back(std::move(new_neighbour_block));

    clusters_to_ix_mapping.emplace_back(-1);
  }
  
  if(restreaming) {
    // Update cluster volume if the streamed node has moved to another cluster during the restreaming phase.  
    if (decision != (*config.stream_nodes_assign)[curr_node_id]) {
        for(PartitionID ix = 0; ix < neighbours.size(); ix++) {
          if(neighbours[ix].first == -1) {
            break;
          }
          blocks[decision].sigma_cluster += neighbours[ix].second;
          blocks[(*config.stream_nodes_assign)[curr_node_id]].sigma_cluster -= neighbours[ix].second;
      }
    }
  } else {
    set_decision(decision, curr_node_id, curr_node_weight, my_thread, streamed_curr_node_edge_count);
  }

  return decision;
}


inline void
vertex_partitioning::set_parent_block(floating_block *parent_block) {
  this->parent_block = parent_block;
}

inline void vertex_partitioning::set_original_problem(
    vertex_partitioning *original_problem) {
  this->original_problem = original_problem;
}

inline floating_block *
vertex_partitioning::get_block_address(PartitionID block_id) const {
  floating_block *block_pointer = &(blocks[block_id]);
  return block_pointer;
}

inline void vertex_partitioning::enable_self_sorting_array() {
  this->use_self_sorting_array = true;
}

inline void vertex_partitioning::increment_graph_edge_count(EdgeWeight edge_count, int restreaming) {
  if(!restreaming) {
    this->current_graph_edge_count += edge_count;
  }
  this->streamed_curr_node_edge_count += edge_count;
}

inline void vertex_partitioning::reset_streamed_edge_count() {
  this->streamed_curr_node_edge_count = 0;
}

inline void vertex_partitioning::increment_curr_node_edge(EdgeWeight node_edge_count) {
  this->streamed_curr_node_edge_count += node_edge_count;
}

inline void vertex_partitioning::increment_cluster_sigma(PartitionID clusterID) {
  blocks[clusterID].update_cluster_sigma(streamed_curr_node_edge_count);
}

inline double vertex_partitioning::calculate_overall_score(std::vector<NodeID> &blocks_weight, std::vector<std::pair<EdgeWeight,EdgeWeight>> &blocks_sigma, LongEdgeID nmbEdges) {
  return 0;
}

#endif /* end of include guard: VERTEX_PARTITIONING_7I4IR31Y */
