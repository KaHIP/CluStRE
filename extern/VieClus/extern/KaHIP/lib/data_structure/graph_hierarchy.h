/******************************************************************************
 * graph_hierarchy.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef GRAPH_HIERACHY_UMHG74CO_2
#define GRAPH_HIERACHY_UMHG74CO_2

#include <stack>

//#include "graph_access.h"
//#include "uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/partial_boundary.h"

namespace KaHIP {

class graph_hierarchy {
public:
        graph_hierarchy( );
        virtual ~graph_hierarchy();

        void push_back(KaHIP::graph_access * G, CoarseMapping * coarse_mapping);
        
        KaHIP::graph_access  * pop_finer_and_project();
        KaHIP::graph_access  * pop_finer_and_project_ns( PartialBoundary & separator );
        KaHIP::graph_access  * get_coarsest();
        CoarseMapping * get_mapping_of_current_finer();
               
        bool isEmpty();
        unsigned int size();
private:
        //private functions
        KaHIP::graph_access * pop_coarsest();

        std::stack<KaHIP::graph_access*>   m_the_graph_hierarchy;
        std::stack<CoarseMapping*>  m_the_mappings;
        std::vector<CoarseMapping*> m_to_delete_mappings;
        std::vector<KaHIP::graph_access*>  m_to_delete_hierachies;
        KaHIP::graph_access  * m_current_coarser_graph;
        KaHIP::graph_access  * m_coarsest_graph;
        CoarseMapping * m_current_coarse_mapping;
};

}

#endif /* end of include guard: GRAPH_HIERACHY_UMHG74CO */
