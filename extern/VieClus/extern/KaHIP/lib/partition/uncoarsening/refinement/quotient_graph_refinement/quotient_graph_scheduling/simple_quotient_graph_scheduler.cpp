/******************************************************************************
 * simple_quotient_graph_scheduler.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/

//#include "random_functions.h"
//#include "simple_quotient_graph_scheduler.h"

#include "extern/KaHIP/lib/tools/random_functions.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/quotient_graph_scheduling/simple_quotient_graph_scheduler.h"

simple_quotient_graph_scheduler::simple_quotient_graph_scheduler(KaHIP::PartitionConfig & config, 
                                                                 QuotientGraphEdges & qgraph_edges,  
                                                                 unsigned int account) {
        unsigned added_edges = 0;
        for( unsigned i = 0; i < (unsigned)ceil(config.bank_account_factor) && added_edges <= account; i++) {
                KaHIP::random_functions::permutate_vector_good_small(qgraph_edges);               
                 for( unsigned i = 0; i < qgraph_edges.size() && added_edges <= account; i++) {
                        m_quotient_graph_edges_pool.push_back(qgraph_edges[i]);
                        added_edges++;
                 }
        }

}

simple_quotient_graph_scheduler::~simple_quotient_graph_scheduler() {
                
}

