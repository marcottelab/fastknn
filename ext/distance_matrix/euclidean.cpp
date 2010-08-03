#include "euclidean.h"


unsigned int unnormalized_manhattan(size_t m, size_t n, size_t k) {
    return (m + n - 2*k);
}

////////////////////////////////////////////////////////////////////////////////
// manhattan
//
// m    defective
// n    drawn
// k    drawn defective
// N, total possible items, does not matter for Manhattan.
//
// Returns a double representing the Euclidean distance between two vectors.
// This is equivalent to the hamming distance.
// Normalizes by the total number of dimensions.
double manhattan(size_t m, size_t n, size_t k, size_t N) {
    return (unnormalized_manhattan(m,n,k)) / (double)(N);
}

////////////////////////////////////////////////////////////////////////////////
// euclidean
//
// m    defective
// n    drawn
// k    drawn defective
// N, total possible items, does not matter for Euclidean.
//
// Returns a double representing the Euclidean distance between two vectors.
// NOTE: Values in the vectors must be 1 or 0.
double euclidean(size_t m, size_t n, size_t k, size_t N) {
    return sqrt(unnormalized_manhattan(m, n, k)) / (double)(N);
}
