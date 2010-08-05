/*
 * File:   similarity.cpp
 * Author: jwoods
 *
 * Created on August 1, 2010, 7:25pm
 */

#include "similarity.h"


////////////////////////////////////////////////////////////////////////////////
// jaccard
//
// m    defective (dimensional length = 1)
// n    drawn (dimensional length = 1)
// k    drawn defective (dimensional length = 1 for both)
// N    total possible items (UNUSED)
//
// Returns a double representing the jaccard distance between the two sets.
// Jaccard distance doesn't care about the total number of genes!
double jaccard(size_t m, size_t n, size_t k, size_t N) {
    return (m + n - 2*k) / (double)(m + n - k);
}


// TODO: Check that this is the correct definition for hellinger.
////////////////////////////////////////////////////////////////////////////////
// hellinger
//
// m    defective (dimensional length = 1)
// n    drawn (dimensional length = 1)
// k    drawn defective (dimensional length = 1 for both)
// N    total possible items (UNUSED)
//
// Returns a double representing the Hellinger distance (1.0 - Sorenson index)
// of two sets.
// Like Jaccard distance, this doesn't care about the total number of genes!
double sorensen(size_t m, size_t n, size_t k, size_t N) {
    return 1.0 - ((2*k) / (double)(m+n));
}
