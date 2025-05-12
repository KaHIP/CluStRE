/*****************************************************************************
 * streamcpi.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/

#include <argtable3.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <regex.h>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <sys/resource.h>
#include <vector>

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

#include "data_structure/ExternalPQ.h"
#include "data_structure/graph_access.h"
#include "graph_io_stream.h"
#include "macros_assertions.h"
#include "parse_parameters.h"
#include "partition/partition_config.h"
#include "quality_metrics.h"
#include "timer.h"
#include "tools/random_functions.h"

#include "partition/onepass_partitioning/leiden.h"
#include "partition/onepass_partitioning/modularity.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"

#include "extclustering/extclustering.h"
#include "extclustering/extclustering_vieclus.h"

#include "FlatBufferWriter.h"
#include "Stream_CPI_Info_generated.h"
#include "cpi/run_length_compression.hpp"

#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "data_structure/compression_vectors/BatchRunLengthCompression.h"

#include "robin_hood.h"
#include <chrono>
// #include "absl/container/flat_hash_set.h"
// #include "absl/container/flat_hash_map.h"

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

// Struct to store captured values
struct CapturedValues {
    std::size_t space_in_bytes;
    std::size_t uncompressed_space_in_bytes;
    double space_in_mib;
    double relative;
};

void MemoryConsumptionSignificantDS(robin_hood::unordered_flat_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> quotient,
                                    std::vector<floating_block> &artificial_blocks,
                                    std::vector<NodeWeight> &stream_blocks_weight,
                                    std::shared_ptr<CompressionDataStructure<PartitionID>> &block_assignments,
                                    std::vector<PartitionID> &stream_nodes_assigned,
                                    std::vector<PartitionID> &clusters_mapping,
                                    std::vector<std::vector<std::pair<PartitionID, EdgeWeight>>> &neighbor_blocks,
                                    LongNodeID rle_length);

void RunningTimeSubModules(double &global_mapping_time,
                           double &buffer_io_time,
                           double &total_time,
                           double &ext_alg_time,
                           double &atf_node_construction,
                           double &scr_vec_push_back,
                           double &qgraph_update,
                           double &node_assignments,
                           double &io_label_prop_time);

void initialize_onepass_partitioner(HeiClus::PartitionConfig &config,
                                    vertex_partitioning *&onepass_partitioner);

void initialize_extclustering(int argn, char **argv, HeiClus::PartitionConfig &config,
                              extclustering *&ext_clusterer);

std::string extractBaseFilename(const std::string &fullPath);

CapturedValues parseCapturedValues(const std::string &output_str);

long getMaxRSS();

std::ostream &cout_redirect();

std::ostringstream redirected_cout;

int main(int argn, char **argv) {
    HeiClus::PartitionConfig config;
    std::string graph_filename;
    std::vector<std::vector<LongNodeID>> *input = NULL;
    EdgeWeight total_edge_cut = 0;
    EdgeWeight qap = 0;
    double balance = 0;
    double total_nodes = 0;

    int active_nodes_amount = 0;

    timer t, processing_t, io_t, ext_t, atf_t, score_t, qgraph_t, nd_assign_t, io_label_prop_t, ls_frac_time_t, ls_time_t;

    double global_mapping_time = 0;
    double buffer_io_time = 0;
    double total_time = 0;
    double ext_alg_time = 0;
    double atf_node_construction = 0;
    double scr_vec_push_back = 0;
    double qgraph_update = 0;
    double node_assignments = 0;
    double io_lp_time = 0;
    double ls_frac_time = 0;
    double ls_time = 0;

    bool is_graph_weighted = false;
    bool suppress_output = false;
    bool recursive = false;
    bool active_nodes_exist = false;

    // Read and initialise Command line Arguments.
    int ret_code = parse_parameters(argn, argv, config, graph_filename, is_graph_weighted,
                                    suppress_output, recursive);

    if (ret_code) {
        return 0;
    }

    config.stream_input = true;

    // container for storing block assignments used by the streaming algorithm
    std::shared_ptr<CompressionDataStructure<PartitionID>> block_assignments;
    std::cout << R"(
     _____ _        _____ _   _____  ______ 
    / ____| |      / ____| | |  __ \|  ____|
   | |    | |_   _| (___ | |_| |__) | |__   
   | |    | | | | |\___ \| __|  _  /|  __|  
   | |____| | |_| |____) | |_| | \ \| |____ 
    \_____|_|\__,_|_____/ \__|_|  \_\______|
                                           
    )" << std::endl;

    ls_frac_time_t.restart();
    processing_t.restart();

    // Read metadata graph information such as number of nodes, number of edges.
    io_t.restart();
    graph_io_stream::readFirstLineStreamClustering(config, graph_filename, total_edge_cut, qap);
    buffer_io_time += io_t.elapsed();

    // Initialise the selected score function metric for the streaming clustering algorithm
    vertex_partitioning *onepass_partitioner = NULL;
    initialize_onepass_partitioner(config, onepass_partitioner);

    // Initialise the selected in-memory clustering algorithm
    extclustering *ext_clusterer = NULL;
    initialize_extclustering(argn, argv, config, ext_clusterer);

    // set up block assignment container based on algorithm configuration
    if (config.rle_length == 0) {
        block_assignments = std::make_shared<RunLengthCompressionVector<PartitionID>>();
    }

    LongNodeID num_lines = 1;

    for (int restreaming = 0; restreaming < config.restream_amount + 1; restreaming++) {
        
        if(restreaming) {
            
            // Initialse the restreaming phase
            io_t.restart();
            graph_io_stream::restreamingFileReset(config, graph_filename);
            buffer_io_time += io_t.elapsed();

            t.restart();

            while (onepass_partitioner->blocks.size() > (*config.stream_blocks_weight).size()) {
                onepass_partitioner->blocks.pop_back();
                config.clusters_to_ix_mapping.pop_back();
                config.neighbor_blocks[0].pop_back();
            }

            for (PartitionID partition = 0; partition < onepass_partitioner->blocks.size(); partition++) {
                onepass_partitioner->blocks[partition].e_weight = 0;
            }

            global_mapping_time += t.elapsed();

            config.next_key[0] = 0;
        }

        for (LongNodeID curr_node = 0; curr_node < config.n_batches; curr_node++) {
            int my_thread = 0;

            // Load 1 line (1 Node) from disk to then serve an the input to the clusterer
            io_t.restart();
            graph_io_stream::loadBufferLinesToBinary(config, input, num_lines, curr_node, restreaming);
            buffer_io_time += io_t.elapsed();

            t.restart();

            // Insert artificial nodes for neighbouring nodes that have already been streamed
            atf_t.restart();
            graph_io_stream::readNodeOnePassClustering(config, curr_node, restreaming, my_thread, input, block_assignments, onepass_partitioner);
            atf_node_construction += atf_t.elapsed();

            // Calculate score and assign node to cluster with highest score gain
            score_t.restart();
            PartitionID block = onepass_partitioner->solve_node_clustering(
                curr_node, 1, restreaming, my_thread, config.neighbor_blocks[my_thread], config.clusters_to_ix_mapping, config, config.previous_assignment, config.kappa, false);
            scr_vec_push_back += score_t.elapsed();

            // Update Quotient graph for a selected mode
            if (restreaming == 0 && (config.mode != LIGHT && config.mode != LIGHT_PLUS)) {
                qgraph_t.restart();
                onepass_partitioner->update_quotient_graph(block, config.neighbor_blocks[my_thread]);
                qgraph_update += qgraph_t.elapsed();
            }

            PartitionID orig_part = -1;

            // Assign the decision cluster to node.
            nd_assign_t.restart();
            if (config.rle_length == -1) {
                if (restreaming) {
                    orig_part = (*config.stream_nodes_assign)[curr_node];
                }
                (*config.stream_nodes_assign)[curr_node] = block;
            }
            else if (config.rle_length == 0) {
                block_assignments->Append(block);
            }
            node_assignments += nd_assign_t.elapsed();

            // Insert Active nodes when needed
            if (restreaming == config.restream_amount // last streaming phase
                && orig_part != block                 // Active Nodes
                && (config.mode == LIGHT_PLUS || config.mode == STRONG)) {
                active_nodes_exist = true;
                for (LongNodeID &neighbour : (*input)[0]) {
                    if (!config.activeNodes_set->contains(neighbour - 1)) {
                        config.activeNodes_set->insert(neighbour - 1);
                    }
                }
            }

            // Delete currently streamed node from memory
            if (!config.ram_stream) {
                delete input;
            }

            // If new block is the decision, add entry to the stream blocks weight.
            config.previous_assignment = block;

            if (config.stream_blocks_weight->size() <= block) {
                config.stream_blocks_weight->emplace_back(0);
            }
            (*config.stream_blocks_weight)[block] += 1;
            if (restreaming) {
                (*config.stream_blocks_weight)[orig_part] -= 1;
            }
            global_mapping_time += t.elapsed();
        }

        // User defined an external Algorithm to be used.
        if (ext_clusterer != NULL && restreaming == 0 && config.mode != LIGHT && config.mode != LIGHT_PLUS) {
            ext_t.restart();
            ext_clusterer->start_clustering(onepass_partitioner->quotient, onepass_partitioner->blocks, *config.stream_blocks_weight, *config.stream_nodes_assign, (*config.stream_blocks_weight).size(), onepass_partitioner->quotient_edge_count, block_assignments);
            ext_alg_time += ext_t.elapsed();
        }
    }

    if(config.activeNodes_set != NULL) {
        active_nodes_amount = config.activeNodes_set->size();
    } else {
        active_nodes_amount = 0;
    }

    // Local Search Phase
    ls_frac_time += ls_frac_time_t.elapsed();

    if (config.mode == STRONG || config.mode == LIGHT_PLUS) {
        
        // Set timer limit
        auto start_time = std::chrono::steady_clock::now();
        std::chrono::seconds TIME_LIMIT(config.ls_time_limit);

        while (active_nodes_exist) {
            ls_time_t.restart();

            robin_hood::unordered_set<LongNodeID> * new_activeNodes_set = new robin_hood::unordered_set<LongNodeID>();

            onepass_partitioner->curr_round_delta_mod = 0;
            t.restart();
            active_nodes_exist = false;

            for (PartitionID partition = 0; partition < onepass_partitioner->blocks.size(); partition++) {
                onepass_partitioner->blocks[partition].e_weight = 0;
            }
            global_mapping_time += t.elapsed();

            config.next_key[0] = 0;
            LongNodeID previous_node = -1;
            bool timeStop = false;

            for (const LongNodeID &curr_node : *config.activeNodes_set) {

                // Efficiently load neighbours of active node using a saved offset vector
                io_label_prop_t.restart();
                graph_io_stream::loadLinesFromStreamWithOffset(config, input, num_lines, curr_node);
                io_lp_time += io_label_prop_t.elapsed();

                t.restart();
                
                // Insert artificial nodes for neighbouring nodes that have already been streamed
                atf_t.restart();
                graph_io_stream::readNodeOnePassClustering(config, curr_node, 1, 0, input, block_assignments, onepass_partitioner);
                atf_node_construction += atf_t.elapsed();
                
                // Calculate score and assign node to cluster with highest score gain
                score_t.restart();
                PartitionID block = onepass_partitioner->solve_node_clustering(
                    curr_node, 1, 1, 0, config.neighbor_blocks[0], config.clusters_to_ix_mapping, config, config.previous_assignment, config.kappa, true);
                scr_vec_push_back += score_t.elapsed();

                PartitionID orig_part = -1;

                // Assign decision cluster to node.
                nd_assign_t.restart();
                if (config.rle_length == -1) {
                    orig_part = (*config.stream_nodes_assign)[curr_node];
                    (*config.stream_nodes_assign)[curr_node] = block;
                }
                else if (config.rle_length == 0) {
                    block_assignments->Append(block);
                }
                node_assignments += nd_assign_t.elapsed();

                // Add active nodes to hashmap
                if (orig_part != block) {
                    active_nodes_exist = true;
                    for (LongNodeID &neighbour : (*input)[0]) {
                        if(!new_activeNodes_set->contains(neighbour - 1)) {
                            new_activeNodes_set->insert(neighbour - 1);
                        }
                    }
                }

                if (!config.ram_stream) {
                    delete input;
                }

                // If new block is the decision, add entry to the stream blocks weight.
                config.previous_assignment = block;

                if (config.stream_blocks_weight->size() <= block) {
                    config.stream_blocks_weight->emplace_back(0);
                }
                (*config.stream_blocks_weight)[block] += 1;
                (*config.stream_blocks_weight)[orig_part] -= 1;

                global_mapping_time += t.elapsed();

                if (previous_node != -1) {
                    config.activeNodes_set->erase(previous_node);
                }
                previous_node = curr_node;

                // Check if time limit has been reached
                auto elapsed_time = std::chrono::steady_clock::now() - start_time;
                if (elapsed_time > TIME_LIMIT) {
                    timeStop = true;
                    break; // Exit both the for-loop and the while-loop by returning from the function
                }
            }

            delete config.activeNodes_set;
            config.activeNodes_set = new_activeNodes_set;
            new_activeNodes_set = nullptr; 

            // NOTE: Here the overall_delta_mod may not be accurate because it is the approximate modularity computed when streaming
            // When using quotient graph refinement, we also do not factor the modularity upgrade from running VieClus
             if (onepass_partitioner->curr_round_delta_mod < (onepass_partitioner->overall_delta_mod * config.cut_off)) {
                break;
            }
            if (timeStop) {
                active_nodes_exist = false;
            }

            ls_time += ls_time_t.elapsed();
            if(ls_time > config.ls_frac_time * ls_frac_time) {
                break;
            }

            // update the overall delta mod after each round
            onepass_partitioner->overall_delta_mod = onepass_partitioner->overall_delta_mod + onepass_partitioner->curr_round_delta_mod;
        }
    }

    buffer_io_time += io_lp_time;

    if (config.activeNodes_set != NULL) {
        delete config.activeNodes_set;
    }

    total_time += processing_t.elapsed();

    // Later to be used if RAM Streaming is implemented.
    if (config.ram_stream) {
        delete input;
        /* delete lines; */
    }

    long overall_max_RSS = getMaxRSS();
    std::string baseFilename = extractBaseFilename(graph_filename);

    CapturedValues capturedValues;
    
    // write the partition to the disc
    std::stringstream filename;

    std::string mode;
    if (config.mode == LIGHT_PLUS) {
        mode = "light_plus";
    }
    else if (config.mode == LIGHT) {
        mode = "light";
    }
    else if (config.mode == EVO) {
        mode = "evo";
    }
    else if (config.mode == STANDARD) {
        mode = "standard";
    }
    else if (config.mode == STRONG) {
        mode = "strong";
    }

    if (!config.filename_output.compare("")) {
        filename << config.output_path << baseFilename << "_" << mode << ".txt";
    } else {
        filename << config.filename_output;
    }

    if (!config.suppress_file_output) {
        if (ext_clusterer != NULL && config.rle_length == 0) {
            ext_clusterer->writeClusteringToDisk(config, filename.str(), block_assignments);
        } else {
            graph_io_stream::writeClusteringStream(config, filename.str(), onepass_partitioner, block_assignments);
        }
    }

    // Evaluate score of clustering if wished.
    if (config.evaluate) {
        if (ext_clusterer != NULL && config.rle_length == 0) {
            block_assignments = graph_io_stream::readClustering(config, filename.str());
        }
        graph_io_stream::streamEvaluateClustering(config, graph_filename, onepass_partitioner, block_assignments);
    }

    // Get Clusters amount
    int clusters_amount = 0;
    for (auto &block_weight : (*config.stream_blocks_weight)) {
        if (block_weight > 0) {
            clusters_amount++;
        }
    }

    for(auto & block : onepass_partitioner->blocks) {
        if(block.sigma_cluster > config.max_cluster_edge_volume) {
            std::cout<<"error"<<std::endl;
        }
    }

    // Update FlatBuffer file
    FlatBufferWriter fb_writer;
    fb_writer.updateCompressionStatistics(capturedValues.space_in_bytes, capturedValues.uncompressed_space_in_bytes,
                                          capturedValues.space_in_mib, capturedValues.relative);
    fb_writer.updateResourceConsumption(buffer_io_time, ext_alg_time, global_mapping_time,
                                        total_time, overall_max_RSS);
    fb_writer.updateClusteringMetrics(config.score, clusters_amount, active_nodes_amount);
    fb_writer.writeClustering(baseFilename, config);

    return 0;
}

void initialize_onepass_partitioner(HeiClus::PartitionConfig &config,
                                    vertex_partitioning *&onepass_partitioner) {
    switch (config.one_pass_algorithm) {
    case ONEPASS_LEIDEN:
        onepass_partitioner =
            new onepass_leiden(0, 0, config.total_nodes,
                               config.parallel_nodes, config.mode, false, config.cpm_gamma);
        break;
    case ONEPASS_MODULARITY:
        onepass_partitioner =
            new onepass_modularity(0, 0, config.total_nodes,
                                   config.parallel_nodes, config.mode, false, config.cpm_gamma);
        break;
    }
}

void initialize_extclustering(int argn, char **argv, HeiClus::PartitionConfig &config,
                              extclustering *&ext_clusterer) {

    switch (config.ext_clustering_algorithm) {
    case EXT_VIECLUS_ALGORITHM:
        ext_clusterer =
            new extclustering_vieclus(config.rle_length, argn, argv, config);
        break;
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

// Redirect cout to the stringstream
std::ostream &cout_redirect() {
    static std::ostream cout_redirector(redirected_cout.rdbuf());
    return cout_redirector;
}

// Function to parse the captured values from the redirected output
CapturedValues parseCapturedValues(const std::string &output_str) {
    CapturedValues values;

    // Extract values using stream extraction
    std::istringstream stream(output_str);

    // Ignore text up to '=' and then extract values
    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.space_in_bytes;

    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.uncompressed_space_in_bytes;

    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.space_in_mib;

    stream.ignore(std::numeric_limits<std::streamsize>::max(), '=');
    stream >> values.relative;

    return values;
}

void MemoryConsumptionSignificantDS(robin_hood::unordered_flat_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> quotient,
                                    std::vector<floating_block> &artificial_blocks,
                                    std::vector<NodeWeight> &stream_blocks_weight,
                                    std::shared_ptr<CompressionDataStructure<PartitionID>> &block_assignments,
                                    std::vector<PartitionID> &stream_nodes_assigned,
                                    std::vector<PartitionID> &clusters_mapping,
                                    std::vector<std::vector<std::pair<PartitionID, EdgeWeight>>> &neighbor_blocks,
                                    LongNodeID rle_length) {

    std::cout << "----------Memory Consumption Significant DS--------------" << std::endl;

    // Get Memory Consumption of all big initialised vectors, including:
    // 1. Quotient Graph, 2. block weights vector, 3. clusters mapping vector, 4. artificial nodes block 5. block assignments / stream_nodes_assignment

    // assigned nodes vector:
    if (rle_length == -1) {
        size_t AssignedNodesMetadataSize = sizeof(stream_nodes_assigned);
        size_t AssignedNodesElementsSize = stream_nodes_assigned.capacity() * sizeof(PartitionID);
        size_t AssignedNodesTotalMemoryConsumption = AssignedNodesMetadataSize + AssignedNodesElementsSize;
        double AssignedNodesTotalMemoryConsumptionKB = AssignedNodesTotalMemoryConsumption / 1024;

        std::cout << "Stream Assigned Nodes memory consumption: " << AssignedNodesTotalMemoryConsumptionKB << " KB" << std::endl;
    }
    
    // block_assignments Compression Vector
    
    /*size_t BlockAssignmentsMetadataSize = sizeof(block_assignments);
    size_t BlockAssignmentsPointersSize = block_assignments.capacity() * sizeof(std::shared_ptr<CompressionDataStructure<PartitionID>>);

    size_t totalControlBlockSize = 0;
    size_t totalObjectSize = 0;

    for (const auto& ptr : block_assignments) {
        if (ptr) {
            totalControlBlockSize += sizeof(*ptr) + sizeof(std::shared_ptr<CompressionDataStructure<PartitionID>>);
            totalObjectSize += sizeof(CompressionDataStructure<PartitionID>);
        }
    }

    size_t BlockAssignmentstotalMemoryConsumption = BlockAssignmentsMetadataSize + BlockAssignmentsPointersSize + totalControlBlockSize + totalObjectSize;
    double BlockAssignmentstotalMemoryConsumptionKB = totalMemoryConsumption / 1024.0;

    std::cout << "Total memory consumption: " << BlockAssignmentstotalMemoryConsumptionKB << " KB" << std::endl;*/

    // Quotient Graph: robin_hood::unordered_flat_map<std::pair<PartitionID, PartitionID>, EdgeWeight, PairHash> quotient

    std::size_t metadataMemoryHash = sizeof(quotient);

    // Memory used for the hash table storage
    std::size_t storageMemoryHash = quotient.size() * sizeof(std::pair<std::pair<PartitionID, PartitionID>, EdgeWeight>);
    std::size_t totalMemHash = metadataMemoryHash + storageMemoryHash;
    double memoryUsageKB = static_cast<double>(totalMemHash) / 1024.0;

    std::cout << "Quotient Graph memory consumption: " << memoryUsageKB << " KB" << std::endl;

    // block weights vector:
    size_t BlockWeightMetadataSize = sizeof(stream_blocks_weight);
    size_t BlockWeightElementsSize = stream_blocks_weight.capacity() * sizeof(NodeWeight);
    size_t BlockWeightTotalMemoryConsumption = BlockWeightMetadataSize + BlockWeightElementsSize;
    double BlockWeightTotalMemoryConsumptionKB = BlockWeightTotalMemoryConsumption / 1024;

    std::cout << "Stream Blocks Weight memory consumption: " << BlockWeightTotalMemoryConsumptionKB << " KB" << std::endl;

    // Artificial Nodes vector
    size_t ANMetadataSize = sizeof(artificial_blocks);
    size_t ANElementsMemorySize = artificial_blocks.capacity() * sizeof(floating_block);
    size_t ANTotalMemoryConsumption = ANMetadataSize + ANElementsMemorySize;
    double ANTotalMemoryConsumptionKB = ANTotalMemoryConsumption / 1024.0;

    std::cout << "Real Clusters Block total memory consumption: " << ANTotalMemoryConsumptionKB << " KB" << std::endl;

    // Clusters Mapping total Memory Consumption - Very Low
    size_t CMMetadataSize = sizeof(clusters_mapping);
    size_t CMElementsMemorySize = clusters_mapping.capacity() * sizeof(PartitionID);
    size_t CMTotalMemoryConsumption = CMMetadataSize + CMElementsMemorySize;
    double CMTotalMemoryConsumptionKB = CMTotalMemoryConsumption / 1024.0;

    std::cout << "Clusters Mapping total memory consumption: " << CMTotalMemoryConsumptionKB << " KB" << std::endl;

    // neighbours nodes
    size_t NNouterVectorMetadataSize = sizeof(neighbor_blocks);
    size_t NNouterVectorPointersSize = neighbor_blocks.capacity() * sizeof(std::vector<std::pair<PartitionID, EdgeWeight>>);

    size_t NNinnerVectorsMetadataSize = 0;
    size_t NNinnerElementsSize = 0;

    for (const auto &inner_vector : neighbor_blocks) {
        NNinnerVectorsMetadataSize += sizeof(inner_vector);
        NNinnerElementsSize += inner_vector.capacity() * sizeof(std::pair<PartitionID, EdgeWeight>);
    }

    size_t NNtotalMemoryConsumption = NNouterVectorMetadataSize + NNouterVectorPointersSize + NNinnerVectorsMetadataSize + NNinnerElementsSize;

    // Convert to KB
    double NNtotalMemoryConsumptionKB = NNtotalMemoryConsumption / 1024.0;

    // Print results
    std::cout << "Artificial Nodes memory consumption: " << NNtotalMemoryConsumptionKB << " KB" << std::endl;

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void RunningTimeSubModules(double &global_mapping_time,
                           double &buffer_io_time,
                           double &total_time,
                           double &ext_alg_time,
                           double &atf_node_construction,
                           double &scr_vec_push_back,
                           double &qgraph_update,
                           double &node_assignments,
                           double &io_label_prop_time) {

    std::cout << "--------------------Running Time SubModules-----------------" << std::endl;

    std::cout << "Streaming Graph Clustering Time no I/O : " << global_mapping_time << std::endl;
    std::cout << "I/O Time : " << buffer_io_time << std::endl;
    std::cout << "Total Time : " << total_time << std::endl;
    std::cout << "External Algorithm Time : " << ext_alg_time << std::endl;
    std::cout << "Time to build artificial nodes : " << atf_node_construction << std::endl;
    std::cout << "Score and DS push back ops Time : " << scr_vec_push_back << std::endl;
    std::cout << "Quotient Graph Update : " << qgraph_update << std::endl;
    std::cout << "Time to Assign nodes to computed Clusters : " << node_assignments << std::endl;
    std::cout << "IO time only in the label propagation phase : " << io_label_prop_time << std::endl;

    std::cout << "--------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}
