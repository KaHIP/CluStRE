/******************************************************************************
 * buffered_map.h 
 * *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 * Marcelo Fonseca Faraj
 *****************************************************************************/

#ifndef HASHMAP_EFRXO4X2
#define HASHMAP_EFRXO4X2

#include <iostream>
#include <vector>
#include <unordered_map>
#include <bitset>
#include <utility>
#include <algorithm>
#include <set>

#include "definitions.h"
#include <cstdint>

class FNV1aHash {
public:
    inline uint32_t operator()(const int& data) const {
        return fnv1a(data);
    }
};


#endif /* end of include guard: HASHMAP_EFRXO4X2 */
//macros  


