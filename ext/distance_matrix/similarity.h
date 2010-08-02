/*
 * File:   similarity.h
 * Author: jwoods
 *
 * Created on August 1, 2010, 7:25pm
 */

#ifndef SIMILARITY_H_
# define SIMILARITY_H_
double jaccard(size_t, size_t, size_t, size_t);
double hellinger(size_t, size_t, size_t, size_t);
double hamming(size_t, size_t, size_t, size_t);
#endif
