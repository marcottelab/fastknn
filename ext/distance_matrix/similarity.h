/*
 * File:   similarity.h
 * Author: jwoods
 *
 * Created on August 1, 2010, 7:25pm
 */

#ifndef SIMILARITY_H_
# define SIMILARITY_H_

typedef unsigned int size_t;

double jaccard(size_t, size_t, size_t, size_t);
double hellinger(size_t, size_t, size_t, size_t);
#endif
