# CluStRE: Streaming Graph Clustering Algorithm with Multi-Stage Refinement
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## What is **CluStRE**? 
**CluStRE** is a novel streaming graph clustering algorithm that balances computational efficiency with high-quality clustering using multi-stage refinement. Unlike traditional in-memory clustering approaches, **CluStRE** processes graphs in a streaming setting, significantly reducing memory overhead while leveraging re-streaming and evolutionary heuristics to improve solution quality. Our method dynamically constructs a quotient graph, enabling modularity-based optimization while efficiently handling large-scale graphs. We introduce multiple configurations of **CluStRE** to provide trade-offs between speed, memory consumption, and clustering quality. Experimental evaluations demonstrate that **CluStRE** improves solution quality by 89.8%, operates 2.6X faster, and uses less than two-thirds of the memory required by the state-of-the-art streaming algorithm on average. Moreover, our strongest mode enhances solution quality by up to 150% on average. With this, **CluStRE** achieves comparable solution quality to in-memory algorithms, i.e. over 96% of the quality of clustering approaches, including Louvain, effectively bridging the gap between streaming and traditional clustering methods.

## Installation Notes

### Requirements

* C++-14 ready compiler (g++ version 10+)
* CMake
* Argtable (http://argtable.sourceforge.net/)

### Building CluStRE

To build the software, run
```shell
./compile.sh
```

Alternatively, you can use the standard CMake build process.

The resulting binary is located in `deploy/clustre`.

## Running CluStRE

To cluster a graph in METIS format using CluStRE, run

```shell
./clustre <graph filename> --one_pass_algorithm=modularity --ext_clustering_algorithm=VieClus --mode=<mode>
```
By default, the algorithm stores the resulting cluster assignments in a file identified by `<graph_name>_<mode>.txt`. More information about the quality of the clustering, such as solution quality, running time, memory usage, etc. is written to a FlatBuffer binary file identified by `<graph_name>_<mode>.bin`. By adding the flag `--output_path=<directory_path>` one can specify where to store the files. If the cluster assignments file is not wanted, it can be disabled with the `--suppress_file_output` flag.

The `--one_pass_algorithm` flag determines which objective function to optimise when streaming. By setting this variable to `modularity`, the algorithm will assign the streamed nodes to a cluster where the modularity gain is the highest. The modularity objective function is built in, but a user can implement other objective functions and use the desired one. 

The `--ext_clustering_algorithm` flag determines which in-memory algorithm to use on the quotient graph once all nodes have been streamed. The `VieClus` evolutionary algorithm is the built-in in-memory graph clustering algorithm, however, as before, a user could integrate and implement other in-memory graph clustering algorithms and use the desired one. By default, the VieClus algorithm runs for an estimated 5 minutes. However, this can be changed by adding the flag `--ext_algorithm_time=<time_in_seconds>`.

By adding in the flag `--evaluate` the solution quality will be calculated and stored in the FlatBuffer file.

The `--mode` flag can be set to various values depending on which mode you wish to select. Each mode clusters the graph differently using different phases. Refer to the following table. 

| mode       | Description                                                         |
|------------|---------------------------------------------------------------------|
| light      | One-pass Streaming                                   |
| light_plus | Streaming + Restream with Local Search |
| evo        | Streaming + Memetic Quotient Graph Clustering |
| strong     | Streaming + Memetic Quotient Graph Clustering + Restream with Local Search |

When using the `light_plus` or `strong` mode, the flags `--cut_off=<cut_off_value>` and `--ls_time_limit=<time in seconds>` can restrict the local search phase. By using the --cut_off flag , where cut_off_value is between [0,1], the local search phase will stop when the current round has computed a modularity gain smaller than cut_off_value * estimated_overall_modularity value. The local search could also be time restricted using the `--ls_time_limit=<time in seconds>` flag. The default value for the --cut_off flag is 0.05 and for the --ls_time_limit the default value lies at 600 seconds.

Furthermore, the number of clusters initialised could be capped by using the flag `--cluster_fraction=<cluster_frac_value>`, where cluster_frac_value is between [0,1]. cluster_frac_value * total number of nodes, specifies the maximum number of clusters to initialise during the algorithm. By default, the maximum number of clusters is set to 5 million.

For a complete list of parameters alongside with descriptions, run:

```shell
./clustre --help
```

Note:
- The program stores the results of the executed command in a [flatbuffer](https://github.com/google/flatbuffers) `.bin` file identified by `<graph_name>_<mode>.bin`.
- To Cluster graphs in CluStRE with 64 bit vertex IDs and 64bit edge IDs, edit the CMakeLists.txt file to change `Line 72: option(64BITVERTEXMODE "64 bit mode" OFF)` to
  `Line 72: option(64BITVERTEXMODE "64 bit mode" ON)`, and `Line 65: option(64BITMODE "64 bit mode" OFF)` to `Line 65: option(64BITMODE "64 bit mode" ON)`, then run `./compile.sh`. By default, 64 bit vertex and edge IDs are enabled. 

## Data References
In our work, we performed experiments with graphs sourced from the following repositories:
- SNAP Dataset: https://snap.stanford.edu/data/
- DIMACS Implementation Challenge: https://sites.cc.gatech.edu/dimacs10/index.shtml
- Laboratory for Web Algorithmics: https://law.di.unimi.it/
- Network Repository: https://networkrepository.com/
- Torch Geometric: https://pytorch-geometric.readthedocs.io/en/2.6.0/modules/datasets.html

For our experiments, we converted these graphs to the METIS format, while removing parallel edges, self-loops, and directions, and assigning unitary weight to all nodes and edges. For a description of the METIS graph format, please have a look at the [KaHiP manual](https://github.com/KaHIP/KaHIP/raw/master/manual/kahip.pdf).
