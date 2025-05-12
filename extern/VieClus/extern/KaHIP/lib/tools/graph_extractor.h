/******************************************************************************
 * graph_extractor.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef GRAPH_EXTRACTOR_PDUTVIEF_2
#define GRAPH_EXTRACTOR_PDUTVIEF_2

//#include "data_structure/graph_access.h"
//#include "definitions.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "lib/definitions.h"

namespace KaHIP {

class graph_extractor {
        public:
                graph_extractor();
                virtual ~graph_extractor();

                void extract_block(KaHIP::graph_access & G, 
                                   KaHIP::graph_access & extracted_block, 
                                   PartitionID block, 
                                   std::vector<NodeID> & mapping);

                void extract_two_blocks(KaHIP::graph_access & G, 
                                        KaHIP::graph_access & extracted_block_lhs, 
                                        KaHIP::graph_access & extracted_block_rhs, 
                                        std::vector<NodeID> & mapping_lhs, 
                                        std::vector<NodeID> & mapping_rhs,
                                        NodeWeight & partition_weight_lhs,
                                        NodeWeight & partition_weight_rhs); 

               void extract_two_blocks_connected(KaHIP::graph_access & G, 
                                                 std::vector<NodeID> lhs_nodes,
                                                 std::vector<NodeID> rhs_nodes,
                                                 PartitionID lhs, 
                                                 PartitionID rhs,
                                                 KaHIP::graph_access & pair,
                                                 std::vector<NodeID> & mapping) ;


};
}


#endif /* end of include guard: GRAPH_EXTRACTOR_PDUTVIEF */
