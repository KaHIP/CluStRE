/******************************************************************************
 * construct_partition.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef CONSTRUCT_PARTITION_E86DQF5S
#define CONSTRUCT_PARTITION_E86DQF5S

//#include "data_structure/graph_access.h"
//#include "parallel_mh/population.h"
//#include "partition_config.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/parallel_mh/population.h"
#include "extern/KaHIP/lib/partition/partition_config.h"

class construct_partition {
public:
        construct_partition();
        virtual ~construct_partition();

        void construct_starting_from_partition( KaHIP::PartitionConfig & config, KaHIP::graph_access & G);
        void createIndividuum( KaHIP::PartitionConfig & config, KaHIP::graph_access & G, 
                               Individuum & ind, bool output); 
};


#endif /* end of include guard: CONSTRUCT_PARTITION_E86DQF5S */
