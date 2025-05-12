/******************************************************************************
 * most_balanced_minimum_cuts.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef MOST_BALANCED_MINIMUM_CUTS_SBD5CS
#define MOST_BALANCED_MINIMUM_CUTS_SBD5CS

//#include "data_structure/graph_access.h"
//#include "partition_config.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/partition/partition_config.h"

class most_balanced_minimum_cuts {
        public:
                most_balanced_minimum_cuts();
                virtual ~most_balanced_minimum_cuts();

                void compute_good_balanced_min_cut( KaHIP::graph_access & residualGraph,
                                                    const KaHIP::PartitionConfig & config,
                                                    NodeWeight & perfect_rhs_weight, 
                                                    std::vector< NodeID > & new_rhs_node ); 

        private:
                void build_internal_scc_graph( KaHIP::graph_access & residualGraph,  
                                               std::vector<int> & components, 
                                               int comp_count, 
                                               KaHIP::graph_access & scc_graph);

                void compute_new_rhs( KaHIP::graph_access & scc_graph, 
                                      const KaHIP::PartitionConfig & config,
                                      std::vector< NodeWeight > & comp_weights,
                                      int comp_of_s,
                                      int comp_of_t,
                                      NodeWeight optimal_rhs_weight,
                                      std::vector<int> & comp_for_rhs); 
};


#endif /* end of include guard: MOST_BALANCED_MINIMUM_CUTS_SBD5CS */
