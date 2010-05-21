#include "euclidean.h"

////////////////////////////////////////////////////////////////////////////////
// manhattan
//
// m    defective
// n    drawn
// k    drawn defective
// N, total possible items, does not matter for Manhattan.
//
// Returns a double representing the Euclidean distance between two vectors.
// NOTE: Values in the vectors must be 1 or 0.
double manhattan(size_t m, size_t n, size_t k, size_t N = 0) {
    return (m - k) + (n - k);
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
double euclidean(size_t m, size_t n, size_t k, size_t N = 0) {
    return sqrt(manhattan(m, n, k));
}
