/******************************************************************************
 * random_functions.cpp 
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *****************************************************************************/

//#include "random_functions.h"

#include "extern/KaHIP/lib/tools/random_functions.h"


namespace KaHIP {
MersenneTwister random_functions::m_mt;
int random_functions::m_seed = 0;

random_functions::random_functions()  {
}

random_functions::~random_functions() {
}
}