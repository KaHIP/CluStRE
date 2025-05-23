/******************************************************************************
 * quality_metrics.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/



#ifndef QUALITY_METRICS_10HC2I5M
#define QUALITY_METRICS_10HC2I5M

#include "data_structure/graph_access.h"
#include "partition_config.h"

namespace HeiClus {

class quality_metrics {
public:
        quality_metrics();
        virtual ~quality_metrics ();
	double balance_full_stream(std::vector<NodeWeight> &stream_blocks_weight);

};

}

#endif /* end of include guard: QUALITY_METRICS_10HC2I5M */
