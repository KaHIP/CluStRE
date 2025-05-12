/******************************************************************************
 * path.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/

#//include "path.h"

#include "extern/KaHIP/lib/partition/coarsening/matching/gpa/path.h"

path::path() : head(UNDEFINED_NODE), tail(UNDEFINED_NODE), length(0), active(false) {
                
}

path::path(const NodeID & v) : head(v), tail(v), length(0), active(true) {
                
}

path::~path() {
                
}

