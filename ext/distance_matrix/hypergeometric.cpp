/* 
 * File:   hypergeometric.cpp
 * Author: jwoods
 *
 * Created on July 15, 2009, 2:46 PM
 */

#include "hypergeometric.h"


double factorial_stirling(const size_t& val) {
    return -val + log(2.5066282746271) + (val+0.5) * log(val);
}

double factorial_stirling(const double& val) {
    return -val + log(2.5066282746271) + (val+0.5) * log(val);
}


double quick_ln_fac(const size_t& v) {
    // Use a Stirling approximation above 1000, otherwise use ln_factorial table.

    // Set up the factorial table if necessary
    static vector<double> fac_table(1001, 0.0);
    static bool firstRun = true;

    if (v <= 1000) {

        // On the first run, calculate the factorial table.
        if (firstRun) {
            // cout << "First run quick_ln_fac" << endl;

            double x = fac_table[0] = fac_table[1] = 0.0;

            for (size_t i = 2; i < 1001; ++i) {
                x += log(i);
                fac_table[i] = x;
            }
            firstRun = false;
        }

#ifdef DEBUG
        assert(fac_table[3] != 0.0);
#endif

        return fac_table[v];
    } else {
        // Stirling formula
        return factorial_stirling(v);
    }
    
}


double calc_rs_prob_cyn2(size_t n, size_t k, size_t m, size_t N) {
    return quick_ln_fac(n)    + quick_ln_fac(N-n)     + quick_ln_fac(m) +
           quick_ln_fac(N-m)  - quick_ln_fac(n-k)     - quick_ln_fac(k)
          -quick_ln_fac(m-k)  - quick_ln_fac(N-n-m+k) - quick_ln_fac(N);
}


////////////////////////////////////////////////////////////////////////////////
// hypergeometric
//
// m    defective
// n    drawn
// k    drawn defective
// N    total possible items
//
// Returns a double representing the integration of the hypergeometric probability.
double hypergeometric(size_t m, size_t n, size_t k, size_t N) {

    // Keep track of hypergeometric probs we've already calculated.
    // Use a hash to keep track
    static unordered_map<string,double> seen;

    ostringstream key;
    key << N << ',' << n << ',' << m << ',' << k;

    unordered_map<string,double>::const_iterator seen_it = seen.find(key.str());
    if (seen_it != seen.end()) {
        // cout << "Using hash value." << endl;
        return seen_it->second;
    }

#ifdef DEBUG
    assert(k <= n);
    assert(k <= m);
    assert(n <= N);
    assert(m <= N);
#endif

    size_t min = std::min(m,n);

    double sump = 0.0;
    for (size_t i = k; i <= min; ++i) {
        double p = calc_rs_prob_cyn2(n, i, m, N);
        sump += exp(p);
    }

    if (sump < 0.0)
        sump = 0.0;
    else if (sump > 1.0)
        sump = 1.0;

    // Store in the hash.
    seen[key.str()] = sump;

    return sump;
}
