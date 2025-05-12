/******************************************************************************
 * quality_metrics.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef QUALITY_METRICS_10HC2I5M_2
#define QUALITY_METRICS_10HC2I5M_2

//#include "data_structure/graph_access.h"
//#include "data_structure/matrix/matrix.h"
//#include "partition_config.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/data_structure/matrix/matrix.h"
#include "extern/KaHIP/lib/partition/partition_config.h"

namespace KaHIP {

class quality_metrics {
public:
        quality_metrics();
        virtual ~quality_metrics ();

        EdgeWeight edge_cut(KaHIP::graph_access & G);
        EdgeWeight edge_cut(KaHIP::graph_access & G, int * partition_map); 
        EdgeWeight edge_cut(KaHIP::graph_access & G, PartitionID lhs, PartitionID rhs);
        EdgeWeight max_communication_volume(KaHIP::graph_access & G);
        EdgeWeight min_communication_volume(KaHIP::graph_access & G);
        EdgeWeight max_communication_volume(KaHIP::graph_access & G, int * partition_map);
        EdgeWeight total_communication_volume(KaHIP::graph_access & G); 
        EdgeWeight objective(const KaHIP::PartitionConfig & config, KaHIP::graph_access & G, int * partition_map);
        EdgeWeight edge_cut_connected(KaHIP::graph_access & G, int * partition_map);
        int boundary_nodes(KaHIP::graph_access & G);
        NodeWeight separator_weight(KaHIP::graph_access& G);
        double balance(KaHIP::graph_access & G);
        double balance_edges(KaHIP::graph_access & G);
        double balance_separator(KaHIP::graph_access & G);

        NodeWeight total_qap(KaHIP::graph_access & C, matrix & D, std::vector< NodeID > & rank_assign);
        NodeWeight total_qap(matrix & C, matrix & D, std::vector< NodeID > & rank_assign);
};
}
#endif /* end of include guard: QUALITY_METRICS_10HC2I5M */
