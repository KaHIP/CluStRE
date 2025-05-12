/******************************************************************************
 * evo_clustering.cpp
 *
 * Source of VieClus -- Vienna Graph Clustering 
 *****************************************************************************/

#include <argtable3.h>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h> 

#include "extern/VieClus/lib/VieClus_FlatBufferWriter.h"


/*#include "algorithms/cycle_search.h"
#include "balance_configuration.h"
#include "data_structure/graph_access.h"
#include "graph_io.h"
#include "macros_assertions.h"
#include "parallel_mh_clustering/parallel_mh_async_clustering.h"
#include "parse_parameters.h"
#include "partition/graph_partitioner.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "random_functions.h"
#include "timer.h"*/

#include "extern/KaHIP/lib/algorithms/cycle_search.h"
#include "extern/KaHIP/app/balance_configuration.h"
#include "extern/KaHIP/lib/data_structure/graph_access.h"
#include "extern/KaHIP/lib/io/graph_io.h"
#include "lib/tools/macros_assertions.h"
#include "extern/KaHIP/lib/parallel_mh_clustering/parallel_mh_async_clustering.h"
#include "extern/VieClus/app/parse_parameters.h"
#include "extern/KaHIP/lib/partition/graph_partitioner.h"
#include "extern/KaHIP/lib/partition/partition_config.h"
#include "extern/KaHIP/lib/tools/quality_metrics.h"
#include "extern/KaHIP/lib/tools/random_functions.h"
#include "lib/tools/timer.h"

long getMaxRSS();

std::string extractBaseFilename(const std::string &fullPath);

int main(int argn, char **argv) {

        MPI_Init(&argn, &argv);    /* starts MPI */

        KaHIP::PartitionConfig partition_config;
        std::string graph_filename;
        bool is_graph_weighted = false;
        bool suppress_output   = false;
        bool recursive         = false;
        
        timer io_t, mapping_t, total_t;
        double io_time = 0;
        double mapping_time = 0;
        double total_time = 0;

        int ret_code = parse_parameters(argn, argv, 
                        partition_config, graph_filename, 
                        is_graph_weighted, suppress_output, 
                        recursive); 

        if(ret_code) {
                return 0;
        }

        total_t.restart();

        KaHIP::graph_access G;

        io_t.restart();
        graph_io::readGraphWeighted(G, graph_filename);
        io_time += io_t.elapsed();
        
        omp_set_num_threads(1);
        
        partition_config.k = 1;

        mapping_t.restart();
        //parallel_mh_async_clustering mh;
        //mh.perform_partitioning(partition_config, G);
        
        KaHIP::PartitionConfig copy = partition_config;
        
        copy.upper_bound_partition                     = G.number_of_nodes() + 1;
        copy.graph_allready_partitioned                = false;
        copy.node_ordering                             = RANDOM_NODEORDERING; //Louvain algorithm uses random node ordering
        copy.lm_number_of_label_propagation_levels = 0; //We want to skip the LPP phase in the performClustering method (only louvain is wished)

        LouvainMethod{}.performClustering(copy, &G, true);

        mapping_time += mapping_t.elapsed();

        int rank, size;
        MPI_Comm communicator = MPI_COMM_WORLD; 
        MPI_Comm_rank( communicator, &rank);
        MPI_Comm_size( communicator, &size);

        if( rank == ROOT ) {
            total_time = total_t.elapsed();
            G.set_partition_count(G.get_partition_count_compute());
            //std::cout << "modularity \t\t\t" << ModularityMetric::computeModularity(G) << std::endl;
            //std::cout<<"Memory Consumption" << overall_max_RSS <<std::endl;
       
            long overall_max_RSS = getMaxRSS();
        
            double score = ModularityMetric::computeModularity(G);
            int clusters_amount = G.get_partition_count_compute();

            EdgeID num_edges = G.number_of_edges() / 2;
            EdgeID num_nodes = G.number_of_nodes();

            std::string baseFilename = extractBaseFilename(graph_filename);

            // write the partition to the disc 
            std::string filename;
            if(!partition_config.suppress_file_output) {
                filename = partition_config.output_path + baseFilename + ".txt";
                graph_io::writePartition(G, filename, overall_max_RSS);
            }

            FlatBufferWriter fb_writer;
            fb_writer.updateResourceConsumption(io_time, mapping_time,total_time, overall_max_RSS);
            fb_writer.updateClusteringMetrics(score, clusters_amount);
            fb_writer.writeClustering(baseFilename, partition_config, num_edges ,num_nodes);
        }

        MPI_Finalize();
}



long getMaxRSS() {
    struct rusage usage;

    if (getrusage(RUSAGE_SELF, &usage) == 0) {
        // The maximum resident set size is in kilobytes
        return usage.ru_maxrss;
    } else {
        std::cerr << "Error getting resource usage information." << std::endl;
        // Return a sentinel value or handle the error in an appropriate way
        return -1;
    }
}

// Function to extract the base filename without path and extension
std::string extractBaseFilename(const std::string &fullPath) {
    size_t lastSlash = fullPath.find_last_of('/');
    size_t lastDot = fullPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
        // Found a slash, extract the substring after the last slash
        return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
    }
    else {
        // No slash found, just extract the substring before the last dot
        return fullPath.substr(0, lastDot);
    }
}
