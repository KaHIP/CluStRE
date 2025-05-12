/******************************************************************************
 * coarsening.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef COARSENING_UU97ZBTR
#define COARSENING_UU97ZBTR

//#include "data_structure/graph_access.h"
//#include "data_structure/graph_hierarchy.h"
//#include "partition_config.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/data_structure/graph_hierarchy.h"
#include "extern/KaHIP/lib/partition/partition_config.h"

class coarsening {
public:
        coarsening ();
        virtual ~coarsening ();

        void perform_coarsening(const KaHIP::PartitionConfig & config, KaHIP::graph_access & G, KaHIP::graph_hierarchy & hierarchy);
};

#endif /* end of include guard: COARSENING_UU97ZBTR */
