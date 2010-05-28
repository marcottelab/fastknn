/* 
 * File:   hypergeometric.h
 * Author: jwoods
 *
 * Created on July 15, 2009, 2:46 PM
 */

#ifndef HYPERGEOMETRIC_H_
# define HYPERGEOMETRIC_H_


#include <boost/unordered_map.hpp>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
using boost::unordered_map;
using std::ostringstream;
using std::string;
using std::exp;
using std::log;
using std::vector;


double factorial_stirling(const size_t&);
double factorial_stirling(const double&);
double quick_ln_fac(const size_t&);
double calc_rs_prob_cyn2(size_t, size_t, size_t, size_t);
double hypergeometric(size_t, size_t, size_t, size_t);


#endif
