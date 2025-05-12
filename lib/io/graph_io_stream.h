/******************************************************************************
 * graph_io_stream.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef GRAPHIOSTREAM_H_
#define GRAPHIOSTREAM_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <unordered_map>
#include <vector>

#include "cpi/run_length_compression.hpp"
#include "data_structure/ExternalPQ.h"
#include "data_structure/buffered_map.h"
#include "data_structure/graph_access.h"
#include "data_structure/single_adj_list.h"
#include "definitions.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"
#include "partition/partition_config.h"
#include "random_functions.h"
#include "timer.h"
#include "data_structure/compression_vectors/CompressionDataStructure.h"
#include "data_structure/compression_vectors/RunLengthCompressionVector.h"
#include "data_structure/compression_vectors/BatchRunLengthCompression.h"
//#include <kagen.h>

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

typedef std::vector<std::string> *LINE_BUFFER;

class graph_io_stream {
public:
  graph_io_stream();

  virtual ~graph_io_stream();

  static void restreamingFileReset(HeiClus::PartitionConfig &partition_config,
                                          std::string graph_filename );

  static void readFirstLineStream(HeiClus::PartitionConfig &partition_config,
                                  std::string graph_filename,
                                  EdgeWeight &total_edge_cut, EdgeWeight &qap);
  
  static void readFirstLineStreamClustering(HeiClus::PartitionConfig &partition_config,
                                  std::string graph_filename,
                                  EdgeWeight &total_edge_cut, EdgeWeight &qap);

  static void
  loadRemainingLinesToBinary(HeiClus::PartitionConfig &partition_config,
                             std::vector<std::vector<LongNodeID>> *&input);

  static void 
  loadLinesFromStreamWithOffset(HeiClus::PartitionConfig &partition_config,
                                std::vector<std::vector<LongNodeID>> *&input, 
                                LongNodeID & num_lines,
                                LongNodeID curr_node);

  static std::vector<std::vector<LongNodeID>> * readInputAsGraphBinary( HeiClus::PartitionConfig &partition_config, 
                                                                        LongNodeID & num_lines, 
                                                                        LongNodeID & curr_node,
                                                                        int & restreaming);

  static bool hasEnding(std::string const &string, std::string const &ending);

  static void
  loadBufferLinesToBinary(HeiClus::PartitionConfig &partition_config,
                          std::vector<std::vector<LongNodeID>> *&input,
                          LongNodeID & num_lines,
                          LongNodeID & curr_node,
                          int & restreaming );

  static std::vector<std::vector<LongNodeID>> *
  loadLinesFromStreamToBinary(HeiClus::PartitionConfig &partition_config,
                              LongNodeID & num_lines,
                              LongNodeID & curr_node, 
                              int & restreaming );

  static void
  generateNeighborhood(HeiClus::PartitionConfig &partition_config,
                       std::vector<NodeID> &generated_neighbors,
                       std::vector<std::vector<LongNodeID>> *&input);

  static void readNodeOnePassClustering(HeiClus::PartitionConfig &config, LongNodeID curr_node,
                              int restreaming, int my_thread,
                              std::vector<std::vector<LongNodeID>> *&input,
                              const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments,
                              vertex_partitioning *onepass_partitioner);

  static void streamEvaluateClustering(HeiClus::PartitionConfig &config,
                                       const std::string &filename,
                                       vertex_partitioning *onepass_partitioner,
                                       const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);

  static void writeClusteringStream(HeiClus::PartitionConfig &config,
                                   const std::string &filename,
                                   vertex_partitioning *onepass_partitioner,
                                   const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments);
                     
  static std::shared_ptr<CompressionDataStructure<PartitionID>> readClustering(HeiClus::PartitionConfig &config,
                                                                                const std::string &filename);
};

inline void
graph_io_stream::readNodeOnePassClustering(HeiClus::PartitionConfig &config, LongNodeID curr_node,
                                 int restreaming, int my_thread,
                                 std::vector<std::vector<LongNodeID>> *&input,
                                 const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments,
                                 vertex_partitioning *onepass_partitioner) {
  
  auto &read_ew = config.read_ew;
  auto &read_nw = config.read_nw;
  LongNodeID target;
  NodeWeight weight;
    
  LongNodeID cursor = (config.ram_stream) ? curr_node : 0;

  // 1. Reset neighbour blocks

  auto &next_key = config.next_key[my_thread];
  auto &neighbor_blocks = config.neighbor_blocks[my_thread];
  auto &clusters_to_ix_mapping = config.clusters_to_ix_mapping;
  
  onepass_partitioner->clear_edgeweight_blocks(neighbor_blocks, next_key,
                                               my_thread);

  for (auto & neighbor_block : neighbor_blocks) {
    if(neighbor_block.first == -1) {
      break;
    }
    clusters_to_ix_mapping[neighbor_block.first] = -1;
    neighbor_block.first = -1;
    neighbor_block.second = 0;
  }
  
  onepass_partitioner->reset_streamed_edge_count();
  next_key = 0;

  std::vector<LongNodeID> &line_numbers = (*input)[cursor];
  LongNodeID col_counter = 0;

  weight = (read_nw) ? line_numbers[col_counter++] : 1;
  PartitionID selecting_factor = (1 + (PartitionID)read_ew);
  config.edges = (line_numbers.size() - col_counter) / selecting_factor;
  
  LongNodeID deg = 0;

  // 2. Find neighbouring blocks and weights to the blocks of the current streamed node

  while (col_counter < line_numbers.size()) {
    target = line_numbers[col_counter++];
    deg++;
    EdgeWeight edge_weight = (read_ew) ? line_numbers[col_counter++] : 1;
    onepass_partitioner->increment_graph_edge_count(edge_weight, restreaming);

    // 2.1 Get the cluster of the neighbouring node.

    PartitionID targetGlobalPar, targetGlobalParTest;
    if (config.rle_length == -1) {
      targetGlobalPar = (*config.stream_nodes_assign)[target - 1];
    } else if (config.rle_length == 0) {
      if ((target - 1) < curr_node) {
          targetGlobalPar = block_assignments->GetValueByIndex(target-1);
      } else {
        targetGlobalPar = INVALID_PARTITION;
      }
    }
    
    // 2.2 Check if neighbouring cluster has already been visited yet, if yes update the neighbouring blocks
    if (targetGlobalPar != INVALID_PARTITION) {
      PartitionID key = clusters_to_ix_mapping[targetGlobalPar];
      if (clusters_to_ix_mapping[targetGlobalPar] == -1) {
        clusters_to_ix_mapping[targetGlobalPar] = next_key;
        auto &new_element = neighbor_blocks[next_key];
        new_element.first = targetGlobalPar;
        new_element.second = edge_weight;
        next_key++;   
      } else {
        neighbor_blocks[key].second += edge_weight;
      }
    }
  }

  if (deg > config.max_degree)
    config.max_degree = deg;

  // 3. Save the edges and their weights to the neighbouring cluster ID'S of the current streamed node
  for (int i = 0; i < neighbor_blocks.size(); i++) {
    auto &element = neighbor_blocks[i];
    if(element.first == -1) {
      break;
    }
    onepass_partitioner->load_edge(element.first,
                                   element.second, my_thread);
  }
  
  config.remaining_stream_nodes--;
}


inline void graph_io_stream::loadRemainingLinesToBinary(
    HeiClus::PartitionConfig &partition_config,
    std::vector<std::vector<LongNodeID>> *&input) {
  if (partition_config.ram_stream) {
    /*input = graph_io_stream::loadLinesFromStreamToBinary(
        partition_config, partition_config.remaining_stream_nodes);*/
  }
}

inline void graph_io_stream::loadLinesFromStreamWithOffset( HeiClus::PartitionConfig & partition_config,
                                                            std::vector<std::vector<LongNodeID>> *&input, 
                                                            LongNodeID & num_lines,
                                                            LongNodeID curr_node) {
    
  input = new std::vector<std::vector<LongNodeID>>(num_lines);

  std::string bin_ending(".bin");
  std::string parhip_ending(".parhip");
  if (hasEnding(partition_config.graph_filename, bin_ending) || hasEnding(partition_config.graph_filename, parhip_ending)) {

    LongNodeID batch_size = num_lines;
    LongNodeID node;
    LongNodeID target;
  
    unsigned long long nodes_in_batch = static_cast<unsigned long long>(batch_size);
    unsigned long long *vertex_offsets = new unsigned long long[nodes_in_batch + 1];
    (*partition_config.stream_in).seekg((sizeof(unsigned long long) * 3) + (sizeof(unsigned long long) * curr_node));
    (*partition_config.stream_in).read((char *) (vertex_offsets), (nodes_in_batch + 1) * sizeof(unsigned long long));

    unsigned long long edge_start_pos = vertex_offsets[0];
    unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
    unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
    unsigned long long *edges = new unsigned long long[num_edges_to_read]; // we also need the next vertex offset
    (*partition_config.stream_in).seekg(edge_start_pos);
    (*partition_config.stream_in).read((char *) (edges), (num_edges_to_read) * sizeof(unsigned long long));

    unsigned long long pos = 0;
    for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
      std::vector<LongNodeID> vec;
      input->push_back(vec);
      node = static_cast<NodeID>(i);
      unsigned long long degree = (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
      for (unsigned long long j = 0; j < degree; j++, pos++) {
        LongNodeID target = static_cast<LongNodeID>(edges[pos]);
        (*input)[node].push_back(target + 1);
      }
    }
  
    delete vertex_offsets;
    delete edges;

  } else {
    
    // 1. Find the nearest base line that is a multiple of INTERVAL
    int baseLine = curr_node - (curr_node % partition_config.offset_interval);
    if (baseLine < 0) baseLine = 0; // just in case

    // 2. Get the index into partialOffsets
    int partialIndex = baseLine / partition_config.offset_interval;
    if (partialIndex >= partition_config.partialOffsets->size()) {
        std::cerr << "Error: partial index out of range.\n";
        return 1;
    }

    std::vector<std::string> *lines;
    lines = new std::vector<std::string>(1);
    buffered_input *ss2 = NULL;
    
    // 3. Seek to that position
    (*partition_config.stream_in).seekg((*partition_config.partialOffsets)[partialIndex]);
    if (!(*partition_config.stream_in).good()) {
        std::cerr << "Error: seekg() failed.\n";
        return 1;
    }

    // 4. Sequentially read lines until we reach targetLine
    int currentLine = baseLine;
    while (currentLine <= curr_node) {
      std::getline((*partition_config.stream_in), (*lines)[0]);
      currentLine++;
    }

    ss2 = new buffered_input(lines);
    ss2->simple_scan_line((*input)[0]);

    if (currentLine == curr_node + 1) {
        // `line` now contains line #targetLine
        //std::cout << "Line #" << curr_node + 1 << ": " << (*lines)[0] << "\n";
    } else {
        std::cerr << "Reached EOF or encountered error before line " << curr_node << "\n";
    }
      
    (*lines)[0].clear();
    delete ss2;
    delete lines;
  }
}

inline void graph_io_stream::loadBufferLinesToBinary( HeiClus::PartitionConfig &partition_config,
                                                      std::vector<std::vector<LongNodeID>> *&input, 
                                                      LongNodeID & num_lines, 
                                                      LongNodeID & curr_node,
                                                      int & restreaming ) {
  if (!partition_config.ram_stream) {
    std::string bin_ending(".bin");
    std::string parhip_ending(".parhip");

    // Check wether graph file is a binary file

    if (hasEnding(partition_config.graph_filename, bin_ending) || hasEnding(partition_config.graph_filename, parhip_ending)) {
      input = graph_io_stream::readInputAsGraphBinary(partition_config,
                                                      num_lines, 
                                                      curr_node,
                                                      restreaming);
    } else {
      input = graph_io_stream::loadLinesFromStreamToBinary( partition_config,
                                                            num_lines,
                                                            curr_node,
                                                            restreaming);
    }
  }
}

inline std::vector<std::vector<LongNodeID>> *
graph_io_stream::loadLinesFromStreamToBinary(HeiClus::PartitionConfig & partition_config,
                                             LongNodeID & num_lines,
                                             LongNodeID & curr_node,
                                             int & restreaming) {
  std::vector<std::vector<LongNodeID>> *input;
  input = new std::vector<std::vector<LongNodeID>>(num_lines);
  std::vector<std::string> *lines;
  lines = new std::vector<std::string>(1);
  LongNodeID node_counter = 0;
  buffered_input *ss2 = NULL;
  while (node_counter < num_lines) {
    //initialise the partial offset vector
    if((curr_node % partition_config.offset_interval == 0 || curr_node == 0) 
      && restreaming == 1 
      && (partition_config.mode == LIGHT_PLUS || partition_config.mode == STRONG)) {
      std::streampos pos = (*(partition_config.stream_in)).tellg();
      partition_config.partialOffsets->push_back(pos);
    }
    
    std::getline(*(partition_config.stream_in), (*lines)[0]);
    if ((*lines)[0][0] == '%') { // a comment in the file
      continue;
    }
    
    ss2 = new buffered_input(lines);
    ss2->simple_scan_line((*input)[node_counter++]);
    (*lines)[0].clear();
    delete ss2;
  }
  delete lines;
  return input;
}

inline std::vector<std::vector<LongNodeID>> *
graph_io_stream::readInputAsGraphBinary(  HeiClus::PartitionConfig &partition_config,
                                          LongNodeID & num_lines, 
                                          LongNodeID & curr_node,
                                          int & restreaming) {
  
  std::vector<std::vector<LongNodeID>> *input;
  input = new std::vector<std::vector<LongNodeID>>();

  LongNodeID batch_size = num_lines;
  LongNodeID node;
  LongNodeID target;
  
  unsigned long long nodes_in_batch = static_cast<unsigned long long>(batch_size);
  unsigned long long *vertex_offsets = new unsigned long long[nodes_in_batch + 1];
  (*partition_config.stream_in).seekg(partition_config.bin_start_pos);
  (*partition_config.stream_in).read((char *) (vertex_offsets), (nodes_in_batch + 1) * sizeof(unsigned long long));
  unsigned long long next_pos = partition_config.bin_start_pos + (nodes_in_batch) * sizeof(unsigned long long);

  unsigned long long edge_start_pos = vertex_offsets[0];
  unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
  unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
  unsigned long long *edges = new unsigned long long[num_edges_to_read]; // we also need the next vertex offset
  (*partition_config.stream_in).seekg(edge_start_pos);
  (*partition_config.stream_in).read((char *) (edges), (num_edges_to_read) * sizeof(unsigned long long));
  partition_config.bin_start_pos = next_pos;

  unsigned long long pos = 0;
  for (unsigned long long i = 0; i < nodes_in_batch; ++i) {
    std::vector<LongNodeID> vec;
    input->push_back(vec);
    node = static_cast<NodeID>(i);
    unsigned long long degree = (vertex_offsets[i + 1] - vertex_offsets[i]) / sizeof(unsigned long long);
    for (unsigned long long j = 0; j < degree; j++, pos++) {
      LongNodeID target = static_cast<LongNodeID>(edges[pos]);
      (*input)[node].push_back(target + 1);
    }
  }
  
  delete vertex_offsets;
  delete edges;

  return input;
}


inline void
graph_io_stream::generateNeighborhood(HeiClus::PartitionConfig &partition_config,
                                      std::vector<NodeID> &generated_neighbors,
                                      std::vector<std::vector<LongNodeID>> *&input) {
    input = new std::vector<std::vector<LongNodeID>>(1);
    for(auto & neighbor: generated_neighbors) {
        //std::cout << neighbor << ", ";
        (*input)[0].push_back(neighbor);
    }
}

#endif /*GRAPHIOSTREAM_H_*/
