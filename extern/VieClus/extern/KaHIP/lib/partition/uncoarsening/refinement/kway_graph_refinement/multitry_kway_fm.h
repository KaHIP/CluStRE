/******************************************************************************
 * multitry_kway_fm.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef MULTITRY_KWAYFM_PVGY97EW
#define MULTITRY_KWAYFM_PVGY97EW

#include <vector>

//#include "definitions.h"
//#include "kway_graph_refinement_commons.h"
//#include "uncoarsening/refinement/refinement.h"

#include "lib/definitions.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/refinement.h"

class multitry_kway_fm {
        public:
                multitry_kway_fm( );
                virtual ~multitry_kway_fm();

                int perform_refinement(KaHIP::PartitionConfig & config, KaHIP::graph_access & G, 
                                       KaHIP::complete_boundary & boundary, unsigned rounds, 
                                       bool init_neighbors, unsigned alpha);

                int perform_refinement_around_parts(KaHIP::PartitionConfig & config, KaHIP::graph_access & G, 
                                                    KaHIP::complete_boundary & boundary, bool init_neighbors, 
                                                    unsigned alpha, 
                                                    PartitionID & lhs, PartitionID & rhs,
                                                    std::unordered_map<PartitionID, PartitionID> & touched_blocks);


        private:
                int start_more_locallized_search(KaHIP::PartitionConfig & config, KaHIP::graph_access & G, 
                                                 KaHIP::complete_boundary & boundary, 
                                                 bool init_neighbors, 
                                                 bool compute_touched_blocks, 
                                                 std::unordered_map<PartitionID, PartitionID> & touched_blocks, 
                                                 std::vector<NodeID> & todolist);

                KaHIP::kway_graph_refinement_commons* commons;
};

#endif /* end of include guard: MULTITRY_KWAYFM_PVGY97EW  */


