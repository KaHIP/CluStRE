/******************************************************************************
 * floating_block.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/


#include "floating_block.h"
#include "partition/onepass_partitioning/vertex_partitioning.h"


floating_block::floating_block(vertex_partitioning* parent_problem, PartitionID my_block_id, NodeID n_threads, PartitionID depth/*=0*/) {
	e_weight = 0;
	sigma_cluster = 0;
}


