/******************************************************************************
 * cycle_refinement.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef CYCLE_REFINEMENT_JEPIS3F0
#define CYCLE_REFINEMENT_JEPIS3F0

//#include "advanced_models.h"
//#include "cycle_definitions.h"
//#include "definitions.h"
//#include "greedy_neg_cycle.h"
//#include "random_functions.h"
//#include "uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
//#include "uncoarsening/refinement/refinement.h"

#include "extern/KaHIP/lib/partition/uncoarsening/refinement/cycle_improvements/advanced_models.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/cycle_improvements/cycle_definitions.h"
#include "lib/definitions.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/cycle_improvements/greedy_neg_cycle.h"
#include "extern/KaHIP/lib/tools/random_functions.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/kway_graph_refinement/kway_graph_refinement_commons.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/refinement.h"

class cycle_refinement : public refinement{
        public:
                cycle_refinement();
                virtual ~cycle_refinement();

                EdgeWeight perform_refinement(KaHIP::PartitionConfig & config, 
                                              KaHIP::graph_access & G, 
                                              KaHIP::complete_boundary & boundary);

        private:
                EdgeWeight greedy_ultra_model(KaHIP::PartitionConfig & partition_config, 
                                              KaHIP::graph_access & G, 
                                              KaHIP::complete_boundary & boundary);

                EdgeWeight greedy_ultra_model_plus(KaHIP::PartitionConfig & partition_config, 
                                                   KaHIP::graph_access & G, 
                                                   KaHIP::complete_boundary & boundary);

                EdgeWeight playfield_algorithm(KaHIP::PartitionConfig & partition_config, 
                                               KaHIP::graph_access & G, 
                                               KaHIP::complete_boundary & boundary);

                advanced_models m_advanced_modelling;
};




#endif /* end of include guard: CYCLE_REFINEMENT_JEPIS3F0 */
