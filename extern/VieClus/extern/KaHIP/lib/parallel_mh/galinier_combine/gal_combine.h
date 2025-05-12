/******************************************************************************
 * gal_combine.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef GAL_COMBINE_XDMU5YB7
#define GAL_COMBINE_XDMU5YB7

//#include "partition_config.h"
//#include "data_structure/graph_access.h"

#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/KaHIP/lib/data_structure/graph_access.h"

class gal_combine {
public:
        gal_combine();
        virtual ~gal_combine();

        void perform_gal_combine( KaHIP::PartitionConfig & config, KaHIP::graph_access & G);
};


#endif /* end of include guard: GAL_COMBINE_XDMU5YB7 */
