/* 
 * File:   euclidean.h
 * Author: jwoods
 *
 * Created on November 13, 2009, 3:18 PM
 */

#ifndef EUCLIDEAN_H_
#define	EUCLIDEAN_H_

#include <cmath>
#include <cstddef>

using std::sqrt;

unsigned int unnormalized_manhattan(size_t, size_t, size_t);
double manhattan(size_t, size_t, size_t, size_t);
double euclidean(size_t, size_t, size_t, size_t);
#endif	/* EUCLIDEAN_H_ */

