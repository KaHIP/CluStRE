/******************************************************************************
 * parse_parameters.h
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Christian Schulz <christian.schulz.phone@gmail.com>
 *****************************************************************************/

#ifndef PARSE_PARAMETERS_GPJMGSM8
#define PARSE_PARAMETERS_GPJMGSM8

#include "configuration.h"
#include <omp.h>
#include <sstream>

int parse_parameters(int argn, char **argv, HeiClus::PartitionConfig & partition_config,
                     std::string &graph_filename, bool &is_graph_weighted,
                     bool &suppress_program_output, bool &recursive) {

// Definition of possible paramters to be passed.

  const char *progname = argv[0];

  // Setup argtable parameters.
  struct arg_lit *help = arg_lit0(NULL, "help", "Print help.");
  struct arg_str *filename =
      arg_strn(NULL, NULL, "FILE", 1, 1, "Path to graph file to partition.");
  struct arg_str *filename_output = arg_str0(
      NULL, "output_filename", NULL,
      "Specify the name of the output file (that contains the partition).");
  struct arg_str *output_path = arg_str0(
      NULL, "output_path", NULL, "Specify the path of the output file(s).");
  struct arg_int *user_seed =
      arg_int0(NULL, "seed", NULL, "Seed to use for the PRNG.");
  /*struct arg_int *k =
      arg_int1(NULL, "k", NULL, "Number of blocks to partition the graph.");*/
  struct arg_int *k =
      arg_int0(NULL, "k", NULL, "Number of blocks to partition the graph.");
  struct arg_dbl *imbalance =
      arg_dbl0(NULL, "imbalance", NULL, "Desired balance. Default: 3 (%).");

  struct arg_lit *suppress_output =
      arg_lit0(NULL, "suppress_output", "(Default: output enabled)");
  struct arg_lit *suppress_file_output =
      arg_lit0(NULL, "suppress_file_output", "(Default: file output enabled)");

  struct arg_dbl *time_limit =
      arg_dbl0(NULL, "time_limit", NULL, "Time limit in s. Default 0s .");

  struct arg_lit *balance_edges = arg_lit0(
      NULL, "balance_edges", "Turn on balancing of edges among blocks.");

  // Stream Partition
  struct arg_int *stream_buffer =
      arg_int0(NULL, "stream_buffer", NULL,
               "Buffer size (number of nodes) for stream partitioning or "
               "mapping: Default 32768.");
  struct arg_dbl *kappa = arg_dbl0(NULL, "kappa", NULL,
                                   "Prioritize previous block assignment by "
                                   "kappa times: Default: Disabled(1).");
  struct arg_dbl *cluster_fraction = arg_dbl0(NULL, "cluster_fraction", NULL,
                                   "percentage of nodes as an Upperlimit Value for the amount of Clusters"
                                   "Clusters Amount Upperlimit: Default: 1.");
  struct arg_int *max_clusters =
      arg_int0(NULL, "max_clusters", NULL,
              "The maximum amount of clusters to be initialised in total.");

  struct arg_int *rle_length =
      arg_int0(NULL, "rle_length", NULL,
               "Default = std::vector. Set to 0 for full RLE or any other "
               "value for vector of RLEs.");
  struct arg_int *uncompressed_runs =
      arg_int0(NULL, "uncompressed_runs", NULL,
               "No. of uncompressed runs in cpi vector. Default = 64.");
  struct arg_str *one_pass_algorithm =
      arg_str0(NULL, "one_pass_algorithm", NULL,
               "One-pass partitioning algorithm "
               "(balanced|hashing|greedy|ldgsimple|ldg|fennel|fennelapprosqrt|"
               "chunking|fracgreedy|fennsimplemap). (Default: fennel)");
  struct arg_str *mode =
      arg_str0(NULL, "mode", NULL,
               "Streaming Algorithm mode"
               "(light|light_plus|evo|strong). (Default: strong)");
    //Streaming Graph Clustering
  struct arg_int *ext_algorithm_time =
      arg_int0(NULL, "ext_algorithm_time", NULL, "Time limit in s for the external Algorithm if specified. Default 300s.");
  struct arg_int *ls_time_limit =
      arg_int0(NULL, "ls_time_limit", NULL, "Time limit in seconds for the Local Search phase. Default 600s.");
  struct arg_dbl *ls_frac_time =
      arg_dbl0(NULL, "ls_frac_time", NULL, "The maximum relative amount local search phase should be used");
  struct arg_int *offset_interval = 
        arg_int0(NULL, "offset_interval", NULL,
               "The offset interval number to store line beginnings of txt file to allow fast random access. Default = 10.");
   struct arg_int *restream_amount = 
        arg_int0(NULL, "restream_amount", NULL,
               "Total amount of Graph restreams.");
  struct arg_str *ext_clustering_algorithm =
            arg_str0(NULL, "ext_clustering_algorithm", NULL,
               "external Clustering Algorithm "
               "(VieClus | leiden). (Default: no external algorithm)");
  struct arg_dbl *resolution_param =
      arg_dbl0(NULL, "resolution_param", NULL,
               "Default = 0.05 "
               "value for the resolution Parameter used in the CPM score function");
  struct arg_dbl *cut_off =
      arg_dbl0(NULL, "cut_off", NULL,
               "Default = 0.001 "
               "Cut off value to break out of the local search");
  struct arg_dbl *strong_cut_off =
      arg_dbl0(NULL, "strong_cut_off", NULL,
               "Default = 0.0005 "
               "Cut off value to break out of the local search");
  struct arg_lit *ram_stream = arg_lit0(
      NULL, "ram_stream", "Stream from RAM instead of HD. (Default: disabled)");
  struct arg_lit *evaluate =
      arg_lit0(NULL, "evaluate", "Run evaluator. (Default: disabled)");
  struct arg_lit *set_part_zero =
      arg_lit0(NULL, "set_part_zero",
               "Set part ID of all nodes to 0. (Default: disabled)");
  struct arg_lit *write_results = arg_lit0(
      NULL, "write_results",
      "Write experimental results to flatbuffer file. (Default: disabled)");
  struct arg_lit *stream_output_progress =
      arg_lit0(NULL, "stream_output_progress",
               "Output global partition after each batch is partitioned. "
               "(Default: disabled)");
  struct arg_int *num_streams_passes =
      arg_int0(NULL, "num_streams_passes", NULL,
               "Number of times input graph should be streamed. (Default: 1).");

  // KaGen Streaming
    struct arg_lit *streaming_graph_generation = arg_lit0(
            NULL, "streaming_graph_generation", "Generate graph. (Default: disabled)");
    struct arg_int *nodes_to_generate =
            arg_int0(NULL, "nodes_to_generate", NULL,
                     "Number of nodes to generate. (Default: 0).");
    struct arg_int *kagen_chunk_count =
            arg_int0(NULL, "kagen_chunk_count", NULL,
                     "Number of chunks for streaming graph generation. (Default: 0).");
    struct arg_int *kagen_d_ba =
            arg_int0(NULL, "kagen_d_ba", NULL,
                     "d: average vertex degree for barabassi-albert model. (Default: 0).");
    struct arg_lit *rgg2d = arg_lit0(
            NULL, "rgg2d", "Generate rgg2d graph. (Default: disabled)");
	struct arg_lit *rgg3d = arg_lit0(NULL, "rgg3d", "Generate rgg3d graph. (Default: disabled)"); 
	struct arg_lit *rhg = arg_lit0(NULL, "rhg", "Generate rhg graph. (Default: disabled)"); 
	struct arg_lit *rdg2d = arg_lit0(NULL, "rdg2d", "Generate rdg2d graph. (Default: disabled)"); 
	struct arg_lit *rdg3d = arg_lit0(NULL, "rdg3d", "Generate rdg3d graph. (Default: disabled)"); 
	struct arg_lit *ba = arg_lit0(NULL, "ba", "Generate ba graph. (Default: disabled)"); 

    struct arg_dbl *kagen_r = arg_dbl0(NULL, "kagen_r", NULL,
                                     "r: radius for KaGen (Default: 0).");

    struct arg_dbl *kagen_d_rhg = arg_dbl0(NULL, "kagen_d_rhg", NULL,
                                     "d_rhg: average degree for rhg generator (Default: 0).");
	
    struct arg_dbl *kagen_gamma = arg_dbl0(NULL, "gamma", NULL,
                                     "gamma: power-law exponent for rhg generator (Default: 0).");

  struct arg_end *end = arg_end(100);

  // Define argtable.
  void *argtable[] = {help,
                      filename,
                      user_seed,
#ifdef MODE_FREIGHT_GRAPHS
                      filename_output,
                      suppress_output,
                      ram_stream,
                      evaluate,
                      set_part_zero,
                      write_results,
                      suppress_file_output,
                      rle_length,
                      uncompressed_runs,
                      kappa,
                      cluster_fraction,
                      max_clusters,
                      one_pass_algorithm,
                      mode,
                      offset_interval,
                      restream_amount,
                      ext_algorithm_time,
                      ext_clustering_algorithm,
                      ls_time_limit,
                      ls_frac_time,
                      resolution_param,
                      cut_off,
                      strong_cut_off,
                      output_path,
                      nodes_to_generate,
                      kagen_chunk_count,
                      kagen_d_ba,
					  kagen_d_rhg,
					  kagen_gamma,
                      kagen_r,
                      rgg2d,
					  rgg3d,
					  rdg2d,
					  rdg3d,
					  rhg,
					  ba,

#endif
                      end};
  // Parse arguments.
  int nerrors = arg_parse(argn, argv, argtable);

  // Catch case that help was requested.
  if (help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable, "\n");
    arg_print_glossary(stdout, argtable, "  %-40s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 1;
  }

  if (nerrors > 0) {
    arg_print_errors(stderr, end, progname);
    printf("Try '%s --help' for more information.\n", progname);
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    return 1;
  }
  //make this optional
  if (k->count > 0) {
    partition_config.k = k->ival[0];
  }

  if (filename->count > 0) {
    graph_filename = filename->sval[0];
    partition_config.graph_filename = filename->sval[0];
  }

  recursive = false;

// Assign the passed paramters values to our Partition Configuration Data Structure

  HeiClus::configuration cfg;
// Set the default values of the Partition Configuration  
  cfg.standard(partition_config);
#ifdef MODE_MULTILEVELMAPPING
  cfg.eco(partition_config);
#else
  cfg.strong(partition_config);
#endif
#ifdef MODE_KAFFPA
  cfg.eco(partition_config);
#else
  cfg.strong(partition_config);
#endif

#ifdef MODE_NODESEP
  cfg.eco_separator(partition_config);
#endif

#ifdef MODE_STREAMMULTISECTION
  cfg.stream_map(partition_config);
  partition_config.stream_buffer_len = 1;
#elif defined MODE_FREIGHT_GRAPHS
  cfg.stream_map(partition_config);
  partition_config.stream_buffer_len = 1;
#endif

#ifdef MODE_STREAMMAP
  cfg.stream_map(partition_config);

  if (!hierarchy_parameter_string->count) {
    std::cout << "Please specify the hierarchy using the "
                 "--hierarchy_parameter_string option."
              << std::endl;
    exit(0);
  }
  if (!distance_parameter_string->count) {
    std::cout << "Please specify the distances using the "
                 "--distance_parameter_string option."
              << std::endl;
    exit(0);
  }
#endif

#ifdef MODE_STREAMPARTITION
  cfg.stream_partition(partition_config);
#endif

// Adjust the values of the Partition configuration if values have been passed.

  if (filename_output->count > 0) {
    partition_config.filename_output = filename_output->sval[0];
  }

  if (output_path->count > 0) {
    partition_config.output_path = output_path->sval[0];
  }

  if (imbalance->count > 0) {
    partition_config.epsilon = imbalance->dval[0];
  }

  if (balance_edges->count > 0) {
    partition_config.balance_edges = true;
  }

  if (suppress_output->count > 0) {
    suppress_program_output = true;
    partition_config.suppress_output = true;
  }

  if (suppress_file_output->count > 0) {
    partition_config.suppress_file_output = true;
  }

  if (time_limit->count > 0) {
    partition_config.time_limit = time_limit->dval[0];
  }

  if (user_seed->count > 0) {
    partition_config.seed = user_seed->ival[0];
  }

  if (imbalance->count > 0) {
    partition_config.imbalance = imbalance->dval[0];
    partition_config.stream_global_epsilon =
        (partition_config.imbalance) / 100.;
  }

  if (stream_buffer->count > 0) {
    partition_config.stream_buffer_len = (LongNodeID)stream_buffer->ival[0];
  }

  if (rle_length->count > 0) {
    partition_config.rle_length = (LongNodeID)rle_length->ival[0];
  }

    if (nodes_to_generate->count > 0) {
        partition_config.nodes_to_generate = (LongNodeID)nodes_to_generate->ival[0];
    }

    if (kagen_chunk_count->count > 0) {
        partition_config.kagen_chunk_count = kagen_chunk_count->ival[0];
    }

    if (kagen_d_ba->count > 0) {
        partition_config.kagen_d_ba = kagen_d_ba->ival[0];
    }

  if (uncompressed_runs->count > 0) {
    partition_config.uncompressed_runs = (LongNodeID)uncompressed_runs->ival[0];
  }

  if (kappa->count > 0) {
    partition_config.kappa = kappa->dval[0];
  }
  
  if (cluster_fraction->count > 0) {
    partition_config.cluster_fraction = cluster_fraction->dval[0];
    if(partition_config.cluster_fraction > 1 || partition_config.cluster_fraction <= 0) {
      std::cout<<"cluster fraction must be between (0,1]"<<std::endl;
      return 1;
    }
  }

    if (kagen_r->count > 0) {
        partition_config.kagen_r = kagen_r->dval[0];
    }

    if (kagen_d_rhg->count > 0) {
        partition_config.kagen_d_rhg = kagen_d_rhg->dval[0];
    }

    if (kagen_gamma->count > 0) {
        partition_config.kagen_gamma = kagen_gamma->dval[0];
    }

  if (one_pass_algorithm->count > 0) {
#ifdef MODE_STREAMMULTISECTION
#elif defined MODE_FREIGHT_GRAPHS
#else
    if (!(full_stream_mode->count > 0)) {
      std::cout << "--one_pass_algorithm requires --full_stream_mode."
                << std::endl;
      exit(0);
    }
#endif
    if (strcmp("fennsimplemap", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_FENNEL_MAPSIMPLE;
    } else if (strcmp("balanced", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_BALANCED;
    } else if (strcmp("hashing", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_HASHING;
    } else if (strcmp("hashcrc", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_HASHING_CRC32;
    } else if (strcmp("greedy", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_GREEDY;
    } else if (strcmp("ldgsimple", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_LDG_SIMPLE;
    } else if (strcmp("ldg", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_LDG;
    } else if (strcmp("fennel", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_FENNEL;
    } else if (strcmp("fennelapprosqrt", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_FENNEL_APPROX_SQRT;
    } else if (strcmp("chunking", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_CHUNKING;
    } else if (strcmp("fracgreedy", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_FRACTIONAL_GREEDY;
    } else if (strcmp("leiden", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_LEIDEN;
    } else if (strcmp("modularity", one_pass_algorithm->sval[0]) == 0) {
      partition_config.one_pass_algorithm = ONEPASS_MODULARITY;
    }
  }

  if(ext_clustering_algorithm->count > 0) {
    if (strcmp("VieClus", ext_clustering_algorithm->sval[0]) == 0) {
      partition_config.ext_clustering_algorithm = EXT_VIECLUS_ALGORITHM;
    } else if (strcmp("leiden", ext_clustering_algorithm->sval[0]) == 0) {
      partition_config.ext_clustering_algorithm = EXT_LEIDEN_ALGORITHM;
    } else {
      partition_config.ext_clustering_algorithm = NO_EXT_ALGORITHM;
    }
  }

  if(offset_interval->count > 0) {
    partition_config.offset_interval = offset_interval->ival[0];
  }

  if (ext_algorithm_time->count > 0) {
    partition_config.ext_algorithm_time = ext_algorithm_time->ival[0];
  }

  if (ls_time_limit->count > 0) {
    partition_config.ls_time_limit = ls_time_limit->ival[0];
  }

  if (ls_frac_time->count > 0) {
    partition_config.ls_frac_time = ls_frac_time->dval[0];
  }

  if (max_clusters->count > 0) {
    partition_config.max_num_clusters = max_clusters->ival[0];
  }
  
  
  if(restream_amount->count > 0) {
    partition_config.restream_amount = restream_amount->ival[0];
  }

  if(mode->count > 0) {
    if (strcmp("light_plus", mode->sval[0]) == 0) {
      partition_config.mode = LIGHT_PLUS;
    } else if (strcmp("light", mode->sval[0]) == 0) {
      partition_config.mode = LIGHT;
      partition_config.restream_amount = 0;
    } else if (strcmp("evo", mode->sval[0]) == 0) {
      partition_config.mode = EVO;
      partition_config.restream_amount = 0;
    } else if (strcmp("standard", mode->sval[0]) == 0) {
      partition_config.mode = STANDARD;
    } else if (strcmp("strong", mode->sval[0]) == 0) {
      partition_config.mode = STRONG;
    } else {
      partition_config.mode = STRONG;
    }
  }

  if(resolution_param->count > 0) {
    partition_config.cpm_gamma = resolution_param->dval[0];
  }

  if(cut_off->count > 0) {
    partition_config.cut_off = cut_off->dval[0];
  }

  if(strong_cut_off->count > 0) {
    partition_config.strong_cut_off = strong_cut_off->dval[0];
  }


  if (ram_stream->count > 0) {
    partition_config.ram_stream = true;
  }

    if (rgg2d->count > 0) {
        partition_config.rgg2d = true;
    }

    if (rgg3d->count > 0) {
        partition_config.rgg3d = true;
    }

    if (rdg2d->count > 0) {
        partition_config.rdg2d = true;
    }

    if (rdg3d->count > 0) {
        partition_config.rdg3d = true;
    }

    if (rhg->count > 0) {
        partition_config.rhg = true;
    }

    if (ba->count > 0) {
        partition_config.ba = true;
    }

    if (streaming_graph_generation->count > 0) {
        partition_config.streaming_graph_generation = true;
    }

  if (evaluate->count > 0) {
    partition_config.evaluate = true;
  }

  if (set_part_zero->count > 0) {
    partition_config.set_part_zero = true;
  }

  if (write_results->count > 0) {
    partition_config.write_results = true;
  }

  if (stream_output_progress->count > 0) {
    partition_config.stream_output_progress = true;
  }

  return 0;
}

#endif /* end of include guard: PARSE_PARAMETERS_GPJMGSM8 */
