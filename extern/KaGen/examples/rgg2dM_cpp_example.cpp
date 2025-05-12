#include <kagen.h>
#include <mpi.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <cmath>

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

int compute_max_chunk_size(double r) {
	double bound = 1/(r*r); 
	double chunk_bound = static_cast<unsigned int>(std::floor(bound)); 

	int exponent = 0;
	while(true) {
		unsigned int candidate = 1 << exponent; 
			if (candidate < chunk_bound) {
				exponent++;
				continue; 
			} else {
				exponent--;
				break;
			}
	}
	return 1 << exponent; 
}

double compute_radius(uint64_t edges, unsigned int n) {
	double denom = static_cast<double>(n)*n*1.57; 
	double r = std::sqrt(edges/denom); 
	return r; 
}

int main(int argc, char *argv[]) {

  if (argc < 4) {
    std::cout
        << "Error! Not enough arguments. <vertices> <edges> <seed>"
        << std::endl;
    return 1;
  }

  long long n_in = std::stoll(argv[1]);
  if (n_in < 0 || n_in > std::numeric_limits<unsigned int>::max()) {
		  throw std::out_of_range("Value is out of range for an unsigned int");
  }
  unsigned int n = static_cast<unsigned int>(n_in);
  uint64_t m = std::stoull(argv[2]);

  //double r = std::stod(argv[3]);
  int seed = std::stoi(argv[3]);

	double r = compute_radius(m, n);
	int chunks = compute_max_chunk_size(r); 
	std::cout << "Radius: " << r << ", Max number of chunks available: " << chunks << std::endl;


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

  kagen::StreamingGenerator streamGenerator(MPI_COMM_WORLD, 0);

  streamGenerator.setupConfig_RGG2D_M(n, m, 0, false);
  streamGenerator.setRandomSeed(seed);
  streamGenerator.setupChunkGeneration(MPI_COMM_WORLD);
  overall_time += t.elapsed();
  t.restart();
  unsigned int edges = 0;
  for (unsigned int i = 0; i < n; i++) {
    std::vector<unsigned int> neighbors; 
        streamGenerator.streamVertex(i + 1, MPI_COMM_WORLD, neighbors);
    overall_time += t.elapsed();
    edges += neighbors.size();
	std::cout << i+1 << ": ";
    for (unsigned int j = 0; j < neighbors.size(); j++) {
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
  std::cout << "Number of edges:                 " << edges << std::endl;
  std::cout << "Estimated edges:                 "
            << streamGenerator.estimate_edges() << std::endl;
}
