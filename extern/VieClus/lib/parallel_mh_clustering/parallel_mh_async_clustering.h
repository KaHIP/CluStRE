/******************************************************************************
 * parallel_mh_async_clustering.h
 *
 * Source of VieClus -- Vienna Graph Clustering 
 *****************************************************************************/

#ifndef PARALLEL_MH_ASYNC_HF106Y0G
#define PARALLEL_MH_ASYNC_HF106Y0G

#include <mpi.h>

/*#include "data_structure/graph_access.h"
#include "partition_config.h"
#include "population_clustering.h"
#include "timer.h"*/

#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/VieClus/lib/parallel_mh_clustering/population_clustering.h"
#include "lib/tools/timer.h"

class parallel_mh_async_clustering {
public:
        parallel_mh_async_clustering();
        parallel_mh_async_clustering(MPI_Comm communicator);
        virtual ~parallel_mh_async_clustering();

        void perform_partitioning(const KaHIP::PartitionConfig & graph_partitioner_config, KaHIP::graph_access & G);
        void initialize(KaHIP::PartitionConfig & graph_partitioner_config, KaHIP::graph_access & G);
        double perform_local_partitioning(KaHIP::PartitionConfig & graph_partitioner_config, KaHIP::graph_access & G);
        double collect_best_partitioning(KaHIP::graph_access & G, const KaHIP::PartitionConfig & config);
        void perform_cycle_partitioning(KaHIP::PartitionConfig & graph_partitioner_config, KaHIP::graph_access & G);

private:
        //misc
        const unsigned MASTER;
        timer    m_t;
        int      m_rank;
        int      m_size;
        double   m_time_limit;
        bool     m_termination;
        unsigned m_rounds;

        //the best cut found so far
        PartitionID* m_best_global_map;
        int          m_best_global_objective;
        int          m_best_cycle_objective;

        //island
        population_clustering* m_island;
        MPI_Comm m_communicator;
};


#endif /* end of include guard: PARALLEL_MH_ASYNC_HF106Y0G */
