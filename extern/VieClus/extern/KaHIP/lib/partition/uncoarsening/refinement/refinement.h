/******************************************************************************
 * refinement.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef REFINEMENT_UJN9IBHM
#define REFINEMENT_UJN9IBHM

//#include "data_structure/graph_access.h"
//#include "partition_config.h"
//#include "quotient_graph_refinement/complete_boundary.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/complete_boundary.h"

class refinement {
public:
        refinement( );
        virtual ~refinement();
        
        virtual EdgeWeight perform_refinement(KaHIP::PartitionConfig & config, 
                                              KaHIP::graph_access & G, 
                                              KaHIP::complete_boundary & boundary) = 0;
};


#endif /* end of include guard: REFINEMENT_UJN9IBHM */
