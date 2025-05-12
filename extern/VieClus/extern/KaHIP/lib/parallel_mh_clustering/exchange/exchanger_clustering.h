/******************************************************************************
 * exchanger_clustering.h 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 *****************************************************************************/

#ifndef EXCHANGER_YPB6QKNL
#define EXCHANGER_YPB6QKNL

#include <mpi.h>

//#include "data_structure/graph_access.h"
//#include "parallel_mh_clustering/population_clustering.h"
//#include "partition_config.h"
//#include "tools/quality_metrics.h"

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/parallel_mh_clustering/population_clustering.h"
#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/KaHIP/lib/tools/quality_metrics.h"

class exchanger_clustering {
public:
        exchanger_clustering( MPI_Comm communicator );
        virtual ~exchanger_clustering();

        void diversify_population_clustering( KaHIP::PartitionConfig & config, KaHIP::graph_access & G, population_clustering & island, bool replace );
        void quick_start( KaHIP::PartitionConfig & config,  KaHIP::graph_access & G, population_clustering & island );
        void push_best( KaHIP::PartitionConfig & config,  KaHIP::graph_access & G, population_clustering & island );
        void recv_incoming( KaHIP::PartitionConfig & config,  KaHIP::graph_access & G, population_clustering & island );

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

        double m_prev_best_objective;
        int m_max_num_pushes;
        int m_cur_num_pushes;

        MPI_Comm m_communicator;

        quality_metrics m_qm;
};



#endif /* end of include guard: EXCHANGER_YPB6QKNL */
