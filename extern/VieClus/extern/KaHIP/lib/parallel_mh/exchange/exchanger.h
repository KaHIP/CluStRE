/******************************************************************************
 * exchanger.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef EXCHANGER_YPB6QKNL
#define EXCHANGER_YPB6QKNL

#include <mpi.h>

//#include "data_structure/graph_access.h"
//#include "parallel_mh/population.h"
//#include "partition_config.h"
//#include "tools/quality_metrics.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/parallel_mh/population.h"
#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/KaHIP/lib/tools/quality_metrics.h"

class exchanger {
public:
        exchanger( MPI_Comm communicator );
        virtual ~exchanger();

        void diversify_population( KaHIP::PartitionConfig & config, KaHIP::graph_access & G, population & island, bool replace );
        void quick_start( KaHIP::PartitionConfig & config,  KaHIP::graph_access & G, population & island );
        void push_best( KaHIP::PartitionConfig & config,  KaHIP::graph_access & G, population & island );
        void recv_incoming( KaHIP::PartitionConfig & config,  KaHIP::graph_access & G, population & island );

private:
        void exchange_individum(const KaHIP::PartitionConfig & config, 
                                KaHIP::graph_access & G, 
                                int & from, 
                                int & rank, 
                                int & to, 
                                Individuum & in, Individuum & out);

        std::vector< int* >          m_partition_map_buffers;
        std::vector< MPI_Request* > m_request_pointers;
        std::vector<bool>            m_allready_send_to;

        int m_prev_best_objective;
        int m_max_num_pushes;
        int m_cur_num_pushes;

        MPI_Comm m_communicator;

        quality_metrics m_qm;
};



#endif /* end of include guard: EXCHANGER_YPB6QKNL */
