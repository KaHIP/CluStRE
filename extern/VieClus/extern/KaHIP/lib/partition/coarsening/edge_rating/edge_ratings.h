/******************************************************************************
 * edge_ratings.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef EDGE_RATING_FUNCTIONS_FUCW7H6Y
#define EDGE_RATING_FUNCTIONS_FUCW7H6Y

//#include "data_structure/graph_access.h"
//#include "partition_config.h"     

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/partition/partition_config.h"

class edge_ratings {
public:
        edge_ratings(const KaHIP::PartitionConfig & partition_config);
        virtual ~edge_ratings();

        void rate(KaHIP::graph_access & G, unsigned level);
        void rate_expansion_star_2(KaHIP::graph_access & G);
        void rate_expansion_star(KaHIP::graph_access & G);
        void rate_expansion_star_2_algdist(KaHIP::graph_access & G);
        void rate_inner_outer(KaHIP::graph_access & G);
        void rate_pseudogeom(KaHIP::graph_access & G); 
        void compute_algdist(KaHIP::graph_access & G, std::vector<float> & dist); 
        void rate_separator_addx(KaHIP::graph_access & G);
        void rate_separator_multx(KaHIP::graph_access & G);
        void rate_separator_max(KaHIP::graph_access & G);
        void rate_separator_log(KaHIP::graph_access & G);
        void rate_separator_r1(KaHIP::graph_access & G);
        void rate_separator_r2(KaHIP::graph_access & G);
        void rate_separator_r3(KaHIP::graph_access & G);
        void rate_separator_r4(KaHIP::graph_access & G);
        void rate_separator_r5(KaHIP::graph_access & G);
        void rate_separator_r6(KaHIP::graph_access & G);
        void rate_separator_r7(KaHIP::graph_access & G);
        void rate_separator_r8(KaHIP::graph_access & G);
        void rate_realweight(KaHIP::graph_access & G);

private:
        const KaHIP::PartitionConfig & partition_config;
};


#endif /* end of include guard: EDGE_RATING_FUNCTIONS_FUCW7H6Y */
