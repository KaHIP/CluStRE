/******************************************************************************
 * label_propagation_refinement.h  
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/


#ifndef LABEL_PROPAGATION_REFINEMENT_R4XW141Y
#define LABEL_PROPAGATION_REFINEMENT_R4XW141Y

//#include "definitions.h"
//#include "../refinement.h"

#include "extern/KaHIP/lib/partition/uncoarsening/refinement/refinement.h"
#include "lib/definitions.h"

class label_propagation_refinement : public refinement {
public:
        label_propagation_refinement();
        virtual ~label_propagation_refinement();

        virtual EdgeWeight perform_refinement(KaHIP::PartitionConfig & config, 
                                              KaHIP::graph_access & G, 
                                              KaHIP::complete_boundary & boundary); 
};


#endif /* end of include guard: LABEL_PROPAGATION_REFINEMENT_R4XW141Y */
