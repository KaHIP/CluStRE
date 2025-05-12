/******************************************************************************
 * graph_io_stream.cpp
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#include "graph_io_stream.h"
#include "cpi/run_length_compression.hpp"
#include "timer.h"
#include <math.h>
#include <sstream>
 
//#include <kagen.h>
//#include <mpi.h>

#define MIN(A, B) (((A) < (B)) ? (A) : (B))
#define MAX(A, B) (((A) > (B)) ? (A) : (B))

graph_io_stream::graph_io_stream() {}

graph_io_stream::~graph_io_stream() {}

bool graph_io_stream::hasEnding(std::string const &string, std::string const &ending) {
    if (string.length() >= ending.length()) {
        return (0 == string.compare(string.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void graph_io_stream::restreamingFileReset(HeiClus::PartitionConfig &partition_config,
                                          std::string graph_filename ) {
  if (partition_config.stream_in != NULL) {
    delete partition_config.stream_in;
  }

  std::string bin_ending(".bin");
  std::string parhip_ending(".parhip");
  if (hasEnding(graph_filename, bin_ending) || hasEnding(graph_filename, parhip_ending)) {
      std::vector<unsigned long long> buffer(3, 0);
      partition_config.stream_in = new std::ifstream(graph_filename.c_str(), std::ios::binary | std::ios::in);;
      if ((*partition_config.stream_in)) {
          (*partition_config.stream_in).read((char *) (&buffer[0]), 3 * sizeof(unsigned long long));
      }

      partition_config.bin_start_pos = 3 * sizeof(unsigned long long);

  } else {
    partition_config.stream_in = new std::ifstream(graph_filename.c_str());
    if (!(*(partition_config.stream_in))) {
      std::cerr << "Error opening " << graph_filename << std::endl;
      exit(1);
    }
    std::vector<std::string> * lines;

    lines = new std::vector<std::string>(1);
    std::getline(*(partition_config.stream_in), (*lines)[0]);

    // skip comments
    while ((*lines)[0][0] == '%') {
      std::getline(*(partition_config.stream_in), (*lines)[0]);
    }

    if (partition_config.partialOffsets == NULL &&
       (partition_config.mode == LIGHT_PLUS || partition_config.mode == STRONG)) {
      partition_config.partialOffsets = new std::vector<std::streampos>();
      partition_config.partialOffsets->reserve((partition_config.total_nodes / partition_config.offset_interval) + 1);
    }
    
    delete lines;
  }

  if( partition_config.activeNodes_set == NULL &&
      (partition_config.mode == LIGHT_PLUS || partition_config.mode == STRONG)) {
    partition_config.activeNodes_set = new robin_hood::unordered_set<LongNodeID>();
  }

}

void graph_io_stream::readFirstLineStreamClustering(HeiClus::PartitionConfig &partition_config,
                                          std::string graph_filename,
                                          EdgeWeight &total_edge_cut,
                                          EdgeWeight &qap) {
  
  if (partition_config.stream_in != NULL) {
      delete partition_config.stream_in;
  }

  std::string bin_ending(".bin");
  std::string parhip_ending(".parhip");
  if (hasEnding(graph_filename, bin_ending) || hasEnding(graph_filename, parhip_ending)) {
      std::vector<unsigned long long> buffer(3, 0);
      partition_config.stream_in = new std::ifstream(graph_filename.c_str(), std::ios::binary | std::ios::in);;
      if ((*partition_config.stream_in)) {
          (*partition_config.stream_in).read((char *) (&buffer[0]), 3 * sizeof(unsigned long long));
      }

      unsigned long long version = buffer[0];
      partition_config.remaining_stream_nodes = static_cast<NodeID>(buffer[1]);
      partition_config.remaining_stream_edges = static_cast<NodeID>(buffer[2]) / 2;

      partition_config.bin_start_pos = 3 * sizeof(unsigned long long);
  } else {
    partition_config.stream_in = new std::ifstream(graph_filename.c_str());
    if (!(*(partition_config.stream_in))) {
      std::cerr << "Error opening " << graph_filename << std::endl;
      exit(1);
    }
    std::vector<std::string> *lines;

    lines = new std::vector<std::string>(1);
    std::getline(*(partition_config.stream_in), (*lines)[0]);

    // skip comments
    while ((*lines)[0][0] == '%') {
      std::getline(*(partition_config.stream_in), (*lines)[0]);
    }

    std::stringstream ss((*lines)[0]);
    ss >> partition_config.remaining_stream_nodes;
    ss >> partition_config.remaining_stream_edges;
    ss >> partition_config.remaining_stream_ew;

    delete lines;
  }

  switch (partition_config.remaining_stream_ew) {
  case 1:
    partition_config.read_ew = true;
    break;
  case 10:
    partition_config.read_nw = true;
    break;
  case 11:
    partition_config.read_ew = true;
    partition_config.read_nw = true;
    break;
  }

  partition_config.total_edges = partition_config.remaining_stream_edges;
  partition_config.total_nodes = partition_config.remaining_stream_nodes;

  if(partition_config.max_num_clusters == -1) {
    partition_config.max_num_clusters = partition_config.total_nodes * partition_config.cluster_fraction;
  }

  // Index value INVALID_PARTITION => node has not been streamed yet.
  if (partition_config.stream_nodes_assign == NULL && partition_config.rle_length == -1) {
       partition_config.stream_nodes_assign = new std::vector<PartitionID>(
               partition_config.remaining_stream_nodes, INVALID_PARTITION);
  }

  if (partition_config.stream_blocks_weight == NULL) {
    partition_config.stream_blocks_weight =
        new std::vector<NodeWeight>(0, 0);
  }

  partition_config.total_stream_nodeweight = 0;
  partition_config.total_stream_nodecounter = 0;
  partition_config.stream_n_nodes = partition_config.remaining_stream_nodes;

  if (partition_config.num_streams_passes >
      1 + partition_config.restream_number) {
    partition_config.stream_total_upperbound = ceil(
        ((100 + 1.5 * partition_config.imbalance) / 100.) *
        (partition_config.remaining_stream_nodes / (double)partition_config.k));
  } else {
    partition_config.stream_total_upperbound = ceil(
        ((100 + partition_config.imbalance) / 100.) *
        (partition_config.remaining_stream_nodes / (double)partition_config.k));
  }

  partition_config.fennel_alpha =
      partition_config.remaining_stream_edges *
      std::pow(partition_config.k, partition_config.fennel_gamma - 1) /
      (std::pow(partition_config.remaining_stream_nodes,
                partition_config.fennel_gamma));

  partition_config.fennel_alpha_gamma =
      partition_config.fennel_alpha * partition_config.fennel_gamma;

  if (partition_config.full_stream_mode && !partition_config.restream_number) {
    partition_config.quotient_nodes = 0;
  } else {
    partition_config.quotient_nodes = partition_config.k;
  }

  total_edge_cut = 0;
  qap = 0;
  if (partition_config.stream_buffer_len ==
      0) { // signal of partial restream standard buffer size
    partition_config.stream_buffer_len = (LongNodeID)ceil(
        partition_config.remaining_stream_nodes / (double)partition_config.k);
  }
  partition_config.nmbNodes = MIN(partition_config.stream_buffer_len,
                                  partition_config.remaining_stream_nodes);
  partition_config.n_batches = ceil(partition_config.remaining_stream_nodes /
                                    (double)partition_config.nmbNodes);
  partition_config.curr_batch = 0;
}

void graph_io_stream::streamEvaluateClustering(HeiClus::PartitionConfig &config,
                                              const std::string &filename,
                                              vertex_partitioning *onepass_partitioner,
                                              const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments) {
  
  LongNodeID nmbNodes;
  LongEdgeID nmbEdges;
  int ew = 0;

  std::vector<NodeID>blocks_weights(0);
  std::vector<std::pair<EdgeWeight, EdgeWeight>>blocks_sigma(0); //blocks_sigma[i].first -> small Sigma, blocks_sigma[i].second -> big Sigma 
  NodeWeight total_nodeweight = 0;
  //EdgeWeight total_edgeweight = 0;
  uint64_t total_edgeweight = 0;

  std::string bin_ending(".bin");
  std::string parhip_ending(".parhip");
  std::ifstream *in;

  // Check if input file is a binary
  if (hasEnding(config.graph_filename, bin_ending) || hasEnding(config.graph_filename, parhip_ending)) {
    std::vector<unsigned long long> buffer(3, 0);
    in = new std::ifstream(config.graph_filename.c_str(), std::ios::binary | std::ios::in);;
    if ((*in)) {
      (*in).read((char *) (&buffer[0]), 3 * sizeof(unsigned long long));
    }
    
    unsigned long long version = buffer[0];
    nmbNodes = static_cast<NodeID>(buffer[1]);
    nmbEdges = static_cast<NodeID>(buffer[2]) / 2;

    NodeID remaining_nodes = nmbNodes;
    NodeID batch_size = 1;
    NodeID node = 0;
    NodeID target;
    unsigned long long start_pos = 3 * sizeof(unsigned long long);

    while (node < remaining_nodes) {
      unsigned long long nodes_in_batch = static_cast<unsigned long long>(batch_size);
      unsigned long long *vertex_offsets = new unsigned long long[nodes_in_batch + 1];
      (*in).seekg(start_pos);
      (*in).read((char *) (vertex_offsets), (nodes_in_batch + 1) * sizeof(unsigned long long));
      unsigned long long next_pos = start_pos + (nodes_in_batch) * sizeof(unsigned long long);

      unsigned long long edge_start_pos = vertex_offsets[0];
      unsigned long long num_reads = vertex_offsets[nodes_in_batch] - vertex_offsets[0];
      unsigned long long num_edges_to_read = num_reads / sizeof(unsigned long long);
      unsigned long long *edges = new unsigned long long[num_edges_to_read]; // we also need the next vertex offset
      (*in).seekg(edge_start_pos);
      (*in).read((char *) (edges), (num_edges_to_read) * sizeof(unsigned long long));
      start_pos = next_pos;

      unsigned long long pos = 0;

      PartitionID partitionIDSource;
      if(config.rle_length == -1) {
        partitionIDSource = (*config.stream_nodes_assign)[node];
      } else if (config.rle_length == 0) {
        partitionIDSource = block_assignments->GetValueByIndex(node);
      }

      while(blocks_weights.size() <= partitionIDSource) {
        blocks_weights.emplace_back(0);
        blocks_sigma.emplace_back(std::make_pair(0, 0));
      }

      NodeWeight weight = 1;

      unsigned long long degree = (vertex_offsets[1] - vertex_offsets[0]) / sizeof(unsigned long long);
      for (unsigned long long j = 0; j < degree; j++, pos++) {
        auto target = static_cast<NodeID>(edges[pos]);
          
        EdgeWeight edge_weight = 1;
        total_edgeweight += edge_weight;
          
        PartitionID partitionIDTarget;
        if(config.rle_length == -1) {
          partitionIDTarget = (*config.stream_nodes_assign)[target];
        } else if (config.rle_length == 0) {
          partitionIDTarget = block_assignments->GetValueByIndex(target);
        }

        if (partitionIDSource == partitionIDTarget) {
          blocks_sigma[partitionIDSource].first++;
        }
        blocks_sigma[partitionIDSource].second++;
      }
      delete vertex_offsets;
      delete edges;
      node++;
    }

  } else {
 
    std::vector<std::vector<LongNodeID>> *input;
    std::vector<std::string> *lines;
    lines = new std::vector<std::string>(1);
    LongNodeID node_counter = 0;
    buffered_input *ss2 = NULL;
    std::string line;
    std::ifstream in(filename.c_str());

    if (!in) {
      std::cerr << "Error opening " << filename << std::endl;
      return 1;
    }
  
    std::getline(in, (*lines)[0]);
    while ((*lines)[0][0] == '%') {
      std::getline(in, (*lines)[0]); // a comment in the file
    }
  
    std::stringstream ss((*lines)[0]);
    ss >> nmbNodes;
    ss >> nmbEdges;
    ss >> ew;

    bool read_ew = false;
    bool read_nw = false;
    if (ew == 1) {
      read_ew = true;
    } else if (ew == 11) {
      read_ew = true;
      read_nw = true;
    } else if (ew == 10) {
      read_nw = true;
    }
    
    LongNodeID target;

    while (std::getline(in, (*lines)[0])) {
      if ((*lines)[0][0] == '%') {
        continue; // a comment in the file
      }
      LongNodeID node = node_counter++;

      PartitionID partitionIDSource;
      if(config.rle_length == -1) {
        partitionIDSource = (*config.stream_nodes_assign)[node];
      } else if (config.rle_length == 0) {
        partitionIDSource = block_assignments->GetValueByIndex(node);
      }
      
      input = new std::vector<std::vector<LongNodeID>>(1);
      ss2 = new buffered_input(lines);
      ss2->simple_scan_line((*input)[0]);
      std::vector<LongNodeID> &line_numbers = (*input)[0];
      LongNodeID col_counter = 0;

      while(blocks_weights.size() <= partitionIDSource) {
        blocks_weights.emplace_back(0);
        blocks_sigma.emplace_back(std::make_pair(0, 0));
      }
      blocks_weights[partitionIDSource]++;

      NodeWeight weight = 1;
      if (read_nw) {
        weight = line_numbers[col_counter++];
        total_nodeweight += weight;
      }
      while (col_counter < line_numbers.size()) {
        target = line_numbers[col_counter++];
        target = target - 1;
        EdgeWeight edge_weight = 1;
        
        if (read_ew) {
          edge_weight = line_numbers[col_counter++];
        }
        
        total_edgeweight += edge_weight;
        PartitionID partitionIDTarget;
        if(config.rle_length == -1) {
          partitionIDTarget = (*config.stream_nodes_assign)[target];
        } else if (config.rle_length == 0) {
          partitionIDTarget = block_assignments->GetValueByIndex(target);
        }

        if (partitionIDSource == partitionIDTarget) {
          blocks_sigma[partitionIDSource].first++;
        }
        blocks_sigma[partitionIDSource].second++;
      }

      (*lines)[0].clear();
      delete ss2;
      delete input;
      if (in.eof()) {
        break;
      }
    }
    delete lines;
  }

  config.score = onepass_partitioner->calculate_overall_score(blocks_weights, blocks_sigma, nmbEdges);
}

void graph_io_stream::writeClusteringStream(HeiClus::PartitionConfig &config,
                                           const std::string &filename,
                                           vertex_partitioning *onepass_partitioner,
                                           const std::shared_ptr<CompressionDataStructure<PartitionID>>& block_assignments ) {
  std::ofstream f(filename.c_str());
  //std::cout << "writing cluster to " << filename << " ... " << std::endl;

  if(config.rle_length == -1) {
    for (int node = 0; node < config.total_nodes; node++) {
        f << (*config.stream_nodes_assign)[node] << "\n";
    }
  } else if(config.rle_length == 0) {
    for (int node = 0; node < config.total_nodes; node++) {
        f << block_assignments->GetValueByIndex(node) << "\n";
    }
  }
  
  // robin_hood::unordered_map<std::pair<PartitionID, PartitionID>, EdgeWeight>
  /*PartitionID source = 0;
  for(PartitionID s_cluster = 0; s_cluster < onepass_partitioner->blocks.size(); s_cluster++) {
      f <<"Cluster "<<s_cluster<<" neighbours : ";
      for(PartitionID t_cluster = 0; t_cluster < onepass_partitioner->blocks.size(); t_cluster++) {
        if(s_cluster > t_cluster) {continue;}
        std::pair<PartitionID, PartitionID> key = std::make_pair(s_cluster, t_cluster);
        f << ", Connection to " << t_cluster << " with weight " << (*onepass_partitioner->quotient)[key];
    }
    f << std::endl;
  }*/

  //std::cout<<"Amount of Quotient Edges : "<<onepass_partitioner->quotient_edge_count<<std::endl;

  f.close();

}

std::shared_ptr<CompressionDataStructure<PartitionID>> graph_io_stream::readClustering(HeiClus::PartitionConfig &config,
                                    const std::string &filename) {
  
  std::string line;
  (*config.stream_blocks_weight).clear();
  std::shared_ptr<CompressionDataStructure<PartitionID>> block_assignments = std::make_shared<RunLengthCompressionVector<PartitionID>>();

  // open file for reading
  std::ifstream in(filename.c_str());
  if (!in) {
    std::cerr << "Error opening file" << filename << std::endl;
  }

  for (NodeID node = 0; node < config.total_nodes; node++) {
    // fetch current line
    std::getline(in, line);
    while (line[0] == '%') { // Comments
      std::getline(in, line);
    }
    PartitionID partition = (PartitionID)atol(line.c_str());
    block_assignments->Append(partition);
    while((*config.stream_blocks_weight).size() <= partition) {
      (*config.stream_blocks_weight).emplace_back(0);  
    }
    (*config.stream_blocks_weight)[partition] += 1;
  }

  in.close();
  return block_assignments;
}

