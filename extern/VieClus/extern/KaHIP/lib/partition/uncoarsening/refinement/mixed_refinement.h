/******************************************************************************
 * mixed_refinement.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef MIXED_REFINEMENT_XJC6COP3
#define MIXED_REFINEMENT_XJC6COP3

//#include "definitions.h"
//#include "refinement.h"

#include "lib/definitions.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/refinement.h"

class mixed_refinement : public refinement {
public:
        mixed_refinement( );
        virtual ~mixed_refinement();

        virtual EdgeWeight perform_refinement(KaHIP::PartitionConfig & config, 
                                              KaHIP::graph_access & G, 
                                              KaHIP::complete_boundary & boundary); 
};


#endif /* end of include guard: MIXED_REFINEMENT_XJC6COP3 */
