/*
 * File:   similarity.h
 * Author: jwoods
 *
 * Created on August 1, 2010, 7:25pm
 */

#ifndef SIMILARITY_H_
# define SIMILARITY_H_

#include <cstddef>
using std::size_t;

double jaccard(size_t, size_t, size_t, size_t);
double sorensen(size_t, size_t, size_t, size_t);
#endif
