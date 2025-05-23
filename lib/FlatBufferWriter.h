/******************************************************************************
 * FlatBufferWriter.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Adil Chhabra <adil.chhabra@informatik.uni-heidelberg.de>
 *****************************************************************************/

#ifndef KAHIP_FLATBUFFERWRITER_H
#define KAHIP_FLATBUFFERWRITER_H

#include <fstream>
#include <iostream>
#include <vector>
#include "Stream_CPI_Info_generated.h"
#include "partition/partition_config.h"

class FlatBufferWriter {
private:
    double buffer_io_time_;
    double global_mapping_time_;
    double total_time_;
    double ext_alg_time_;
    long maxRSS_;

    EdgeWeight total_edge_cut_;
    double balance_;

    double score_;
    int clusters_amount_;
    int active_nodes_amount_;

    std::size_t space_in_bytes_;
    std::size_t uncompressed_space_in_bytes_;
    double space_in_mib_;
    double relative_;

public:
    FlatBufferWriter()
            : buffer_io_time_(0.0), global_mapping_time_(0.0),
              total_time_(0.0), maxRSS_(0),
              ext_alg_time_(0.0),
              total_edge_cut_(0), balance_(0.0),
              space_in_bytes_(0), uncompressed_space_in_bytes_(0),
              space_in_mib_(0.0), relative_(0.0),
              active_nodes_amount_(0) {}

    void updateResourceConsumption(double &buffer_io_time,
                                   double &ext_alg_time,
                                   double &global_mapping_time, double &total_time, long &maxRSS) {
        buffer_io_time_ = buffer_io_time;
        global_mapping_time_ = global_mapping_time;
        total_time_ = total_time;
        maxRSS_ = maxRSS;
        ext_alg_time_ = ext_alg_time;
    }

    void updatePartitionMetrics(EdgeWeight &total_edge_cut, double &balance) {
        total_edge_cut_ = total_edge_cut;
        balance_ = balance;
    }

    void updateClusteringMetrics(double &score, int & clusters_amount, int & active_nodes_amount) {
        score_ = score;
        clusters_amount_ = clusters_amount;
        active_nodes_amount_ = active_nodes_amount;
    }

    void updateCompressionStatistics(std::size_t &space_in_bytes, std::size_t &uncompressed_space_in_bytes,
                                     double &space_in_mib, double &relative) {
        space_in_bytes_ = space_in_bytes;
        uncompressed_space_in_bytes_ = uncompressed_space_in_bytes;
        space_in_mib_ = space_in_mib;
        relative_ = relative;
    }

    // Function to extract the base filename without path and extension
    static std::string extractBaseFilename(const std::string &fullPath) {
        size_t lastSlash = fullPath.find_last_of('/');
        size_t lastDot = fullPath.find_last_of('.');

        if (lastSlash != std::string::npos) {
            // Found a slash, extract the substring after the last slash
            return fullPath.substr(lastSlash + 1, lastDot - lastSlash - 1);
        } else {
            // No slash found, just extract the substring before the last dot
            return fullPath.substr(0, lastDot);
        }
    }

    void write(const std::string &baseFilename, const HeiClus::PartitionConfig &config) const {
        // output some information about the partition that we have computed
        flatbuffers::FlatBufferBuilder builder(1024);        
        auto filenameOffset = builder.CreateString(baseFilename);
        StreamCPIInfo::GraphMetadataBuilder metadata_builder(builder);
        metadata_builder.add_filename(filenameOffset);
        metadata_builder.add_num_nodes(config.total_nodes);
        metadata_builder.add_num_edges(config.total_edges);
        metadata_builder.add_max_degree(config.max_degree);
        auto metadata = metadata_builder.Finish();
        //builder.Finish(metadata);
        std::cout << "Graph: " << baseFilename << std::endl;
        std::cout << "Nodes: " << config.total_nodes << std::endl;
        std::cout << "Edges: " << config.total_edges << std::endl;
        std::cout << "Max Deg: " << config.max_degree << std::endl;

        StreamCPIInfo::PartitionConfigurationBuilder config_builder(builder);
        config_builder.add_k(config.k);
        config_builder.add_seed(config.seed);
        config_builder.add_input_balance(config.imbalance);
        if (config.set_part_zero) {
            config_builder.add_set_part_zero(1);
        } else {
            config_builder.add_set_part_zero(0);
        }
        config_builder.add_rle_length(config.rle_length);
        config_builder.add_kappa(config.kappa);
        auto configdata = config_builder.Finish();
        //builder.Finish(configdata);
        
        std::cout << "Blocks (k): " << config.k << std::endl;
        std::cout << "Seed: " << config.seed << std::endl;
        std::cout << "Imbalance: " << config.imbalance << std::endl;
        if (config.set_part_zero) {
            std::cout << "Part ID: 0" << std::endl;
        } else {
            std::cout << "Part ID: Fennel" << std::endl;
        }
        if (config.rle_length > 0) {
            std::cout << "RLE Length: " << config.rle_length << std::endl;
        } else if (config.rle_length == -2) {
            std::cout << "External Memory Priority Queue" << std::endl;
        } else if (config.rle_length == 0) {
            std::cout << "RLE Length: n" << std::endl;
        } else {
            std::cout << "Using std::vector" << std::endl;
        }
        std::cout << "Kappa: " << config.kappa << std::endl;

        StreamCPIInfo::RunTimeBuilder runtime_builder(builder);
        runtime_builder.add_io_time(buffer_io_time_);
        runtime_builder.add_mapping_time(global_mapping_time_);
        runtime_builder.add_total_time(total_time_);
        runtime_builder.add_ext_alg_time(ext_alg_time_);
        auto runtimedata = runtime_builder.Finish();
        builder.Finish(runtimedata);
        std::cout << "IO Time: " << buffer_io_time_ << std::endl;
        std::cout << "Mapping Time: " << global_mapping_time_ << std::endl;
        std::cout << "Total Time: " << total_time_ << std::endl;
        std::cout << "External Algorithm Time: " << ext_alg_time_ << std::endl;

        auto partition_metrics =
                StreamCPIInfo::CreatePartitionMetrics(builder, total_edge_cut_, balance_);
        std::cout << "Edge cut: " << total_edge_cut_ << std::endl;
        std::cout << "Balance: " << balance_ << std::endl;

        auto clustering_metrics =
                StreamCPIInfo::CreateClusteringMetrics(builder, score_, clusters_amount_, active_nodes_amount_);
        std::cout << "Score: " << score_ << std::endl;
        std::cout << "Clustering Amount " << clusters_amount_ << std::endl;
        std::cout << "Active Nodes Amount " << active_nodes_amount_ << std::endl;


        // Create MemoryConsumption
        auto memory_consumption = StreamCPIInfo::CreateMemoryConsumption(
                builder, uncompressed_space_in_bytes_,
                space_in_bytes_, space_in_mib_,
                relative_, maxRSS_);
        if (config.rle_length == -1) {
            std::cout << "Compressed space_in_bytes: " << space_in_bytes_
                      << std::endl;
            std::cout << "Uncompressed_space_in_bytes: "
                      << uncompressed_space_in_bytes_ << std::endl;
            std::cout << "Compressed space_in_mib: " << space_in_mib_
                      << std::endl;
            std::cout << "Relative: " << relative_ << std::endl;
        }
        if (maxRSS_ != -1) {
            std::cout << "Maximum Resident Set Size (KB): " << maxRSS_
                      << std::endl;
        }

        // Create Output
        StreamCPIInfo::PartitionBuilder partition_builder(builder);
        partition_builder.add_graph_metadata(metadata);
        partition_builder.add_partition_configuration(configdata);
        partition_builder.add_runtime(runtimedata);
        partition_builder.add_memory_consumption(memory_consumption);
        partition_builder.add_partition_metrics(partition_metrics);
        auto partition = partition_builder.Finish();
        //builder.Finish(partition);

        //Step 4: Write to File
        if(config.write_results) {
            const uint8_t *bufferPointer = builder.GetBufferPointer();
            int bufferSize = builder.GetSize();
            
            std::string outputFileNameStream;
            outputFileNameStream = config.output_path + baseFilename + "_clustering.bin";
            const char *outputFileName = outputFileNameStream.c_str();
            FILE *file = fopen(outputFileName, "wb");
            fwrite(bufferPointer, 1, bufferSize, file);
            fclose(file);
        }
    }

    void writeClustering(const std::string &baseFilename, const HeiClus::PartitionConfig &config) const {
        // output some information about the partition that we have computed
        flatbuffers::FlatBufferBuilder builder(1024);        

        // Create GraphMetadata
        auto filenameOffset = builder.CreateString(baseFilename);
        StreamCPIInfo::GraphMetadataBuilder metadata_builder(builder);
        metadata_builder.add_filename(filenameOffset);
        metadata_builder.add_num_nodes(config.total_nodes);
        metadata_builder.add_num_edges(config.total_edges);
        metadata_builder.add_max_degree(config.max_degree);
        auto metadata = metadata_builder.Finish();
        std::cout << "Graph: " << baseFilename << std::endl;
        std::cout << "Nodes: " << config.total_nodes << std::endl;
        std::cout << "Edges: " << config.total_edges << std::endl;
        std::cout << "Max Deg: " << config.max_degree << std::endl;

        // Create RunTime
        StreamCPIInfo::RunTimeBuilder runtime_builder(builder);
        runtime_builder.add_io_time(buffer_io_time_);
        runtime_builder.add_mapping_time(global_mapping_time_);
        runtime_builder.add_total_time(total_time_);
        runtime_builder.add_ext_alg_time(ext_alg_time_);
        auto runtimedata = runtime_builder.Finish();
        std::cout << "IO Time: " << buffer_io_time_ << std::endl;
        std::cout << "Mapping Time: " << global_mapping_time_ << std::endl;
        std::cout << "Total Time: " << total_time_ << std::endl;
        std::cout << "External Algorithm Time: " << ext_alg_time_ << std::endl;

        // Create ClusteringMetrics
        auto clustering_metrics = StreamCPIInfo::CreateClusteringMetrics(builder, score_, clusters_amount_, active_nodes_amount_);
        std::cout << "Score: " << score_ << std::endl;
        std::cout << "Clusters Amount: " << clusters_amount_ << std::endl;
        std::cout << "Active Nodes Amount: " << active_nodes_amount_ << std::endl;


        // Create MemoryConsumption
        auto memory_consumption = StreamCPIInfo::CreateMemoryConsumption(
            builder, uncompressed_space_in_bytes_,
            space_in_bytes_, space_in_mib_,
            relative_, maxRSS_);
        if (maxRSS_ != -1) {
            std::cout << "Maximum Resident Set Size (KB): " << maxRSS_ << std::endl;
        }

        StreamCPIInfo::PartitionConfigurationBuilder config_builder(builder);
        config_builder.add_k(config.k);
        config_builder.add_seed(config.seed);
        config_builder.add_input_balance(config.imbalance);
        if (config.set_part_zero) {
            config_builder.add_set_part_zero(1);
        } else {
            config_builder.add_set_part_zero(0);
        }
        config_builder.add_rle_length(config.rle_length);
        config_builder.add_kappa(config.kappa);
        config_builder.add_cluster_frac(config.cluster_fraction);
        auto configdata = config_builder.Finish();

        // Create the final Partition object
        StreamCPIInfo::PartitionBuilder partition_builder(builder);
        partition_builder.add_graph_metadata(metadata);
        partition_builder.add_runtime(runtimedata);
        partition_builder.add_memory_consumption(memory_consumption);
        partition_builder.add_partition_configuration(configdata);
        partition_builder.add_clustering_metrics(clustering_metrics);
        auto partition = partition_builder.Finish();

        // Now, call Finish only once with the final partition object
        builder.Finish(partition);

        // Step 4: Write to File
        const uint8_t *bufferPointer = builder.GetBufferPointer();
        int bufferSize = builder.GetSize();

        std::string outputFileNameStream;

        std::string rle;
        if(config.rle_length == -1) {
            rle = "-1";
        } else {
            rle = std::to_string(config.rle_length);
        }

        std::string mode;
        if(config.mode == LIGHT_PLUS) {
            mode = "light_plus";
            std::cout << "Mode: Light+" << std::endl;
        } else if(config.mode == LIGHT) {
            mode = "light";
            std::cout << "Mode: Light" << std::endl;
        } else if(config.mode == EVO) {
            mode = "evo";
            std::cout << "Mode: Evo" << std::endl;
        } else if(config.mode == STANDARD) {
            mode = "standard";
            std::cout << "Mode: Standard" << std::endl;
        } else if(config.mode == STRONG) {
            mode = "strong";
            std::cout << "Mode: Strong" << std::endl;
        }

        outputFileNameStream = config.output_path + baseFilename + "_" + mode + ".bin";
        const char *outputFileName = outputFileNameStream.c_str();
        FILE *file = fopen(outputFileName, "wb");
        fwrite(bufferPointer, 1, bufferSize, file);
        fclose(file);
    }
};

#endif //KAHIP_FLATBUFFERWRITER_H
