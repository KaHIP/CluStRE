/******************************************************************************
 * wcycle_partitioner.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef WCYCLE_PARTITIONER_EPNDQMK
#define WCYCLE_PARTITIONER_EPNDQMK

/*#include "coarsening/coarsening.h"
#include "coarsening/stop_rules/stop_rules.h"
#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "uncoarsening/refinement/refinement.h"*/

#include "extern/KaHIP/lib/partition/coarsening/coarsening.h"
#include "extern/KaHIP/lib/partition/coarsening/stop_rules/stop_rules.h"
#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/refinement.h"

class wcycle_partitioner {
        public:
                wcycle_partitioner( ) : m_level(0) {};
                virtual ~wcycle_partitioner() {};
                int perform_partitioning( const KaHIP::PartitionConfig & config, 
                                          KaHIP::graph_access & G); 

        private:
                int perform_partitioning_recursive( KaHIP::PartitionConfig & partition_config, 
                                                    KaHIP::graph_access & G, 
                                                    KaHIP::complete_boundary ** c_boundary); 

                unsigned   m_level;
                unsigned   m_deepest_level;
                stop_rule* m_coarsening_stop_rule;

                std::unordered_map<unsigned, bool> m_have_been_level_down;
};

#endif /* end of include guard: WCYCLE_PARTITIONER_EPNDQMK */
