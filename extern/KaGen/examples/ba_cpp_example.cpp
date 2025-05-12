#include <kagen.h>
#include <mpi.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/resource.h>

#include "../timer.h"

long getMaxRSS() {
  struct rusage usage;

  if (getrusage(RUSAGE_SELF, &usage) == 0) {
    return usage.ru_maxrss;
  } else {
    std::cout << "Error getting resource usage information." << std::endl;
    return -1; // maybe handle error different
  }
}

int main(int argc, char *argv[]) {

  if (argc < 5) {
    std::cout
        << "Error! Not enough arguments. <vertices> <chunks> <degree> <seed>"
        << std::endl;
    return 1;
  }

  uint64_t n = std::stoll(argv[1]);
  uint64_t chunks = std::stoll(argv[2]);
  unsigned int d = std::stoi(argv[3]);
  unsigned int seed = std::stoi(argv[4]);

  MPI_Init(&argc, &argv);

  timer t;
  double overall_time = 0;
  /*  std::ofstream out("outStream");

    if (!out.is_open()) {
      std::cout << "Error opening the output file." << std::endl;
      return 1;
    } */

  t.restart();

  kagen::KaGen generator(MPI_COMM_WORLD);

  kagen::StreamingGenerator streamGenerator(MPI_COMM_WORLD, chunks);

  streamGenerator.setupConfig_BA(n, 0, d, false, false);
  streamGenerator.setRandomSeed(seed);
  streamGenerator.setupChunkGeneration(MPI_COMM_WORLD);
  overall_time += t.elapsed();
  t.restart();
  uint64_t edges = 0;
  for (uint64_t i = 0; i < n; i++) {
    std::vector<uint64_t> neighbors;
        streamGenerator.streamVertex(i + 1, MPI_COMM_WORLD, neighbors);
    overall_time += t.elapsed();
    edges += neighbors.size();
	std::cout << i + 1 << ": ";
    for (uint64_t j = 0; j < neighbors.size(); j++) {
		std::cout << neighbors[j] << " ";
    }
	std::cout << std::endl;
    t.restart();
    //  std::cout << neighbors.first + 1 << ": ";
    //  for (int j = 0; j < neighbors.second.size(); j++) {
    //    std::c; out << neighbors.second[j] + 1 << " ";
    //  }
    //  std::cout << std::endl;
  }

  // Generate a RGG2D graph with 16 nodes and 32 edges
  // generator.GenerateRGG2D_Streaming(20, 0.5, false, 4);
  /*
  const kagen::Edgelist&    edges = graph.edges;
  const kagen::VertexRange& range = graph.vertex_range;
  */
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /*
  {
      std::stringstream ss;
      ss << "On PE " << rank << " [" << range.first << ", " << range.second <<
  "): "; for (const auto& [u, v]: edges) { ss << u << "->" << v << " ";
      }
      ss << "\n";
      std::cout << ss.str();
  }
  */

  MPI_Finalize();

  std::cout << "Overall time for RHG generation: " << overall_time << std::endl;
  long maxRSS = getMaxRSS();
  std::cout << "Maximum Resident Set Size (KB):  " << maxRSS << std::endl;
  std::cout << "Number of edges:               " << edges << std::endl;
  std::cout << "Estimated edges:               "
            << streamGenerator.estimate_edges() << std::endl;
}
