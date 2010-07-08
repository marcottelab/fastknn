#include "classifier.h"


NaiveBayes::NaiveBayes(const DistanceMatrix* const rhs, size_t k_, float max_distance_)
        : Classifier(rhs), k(k_), max_distance(max_distance_) { }

// Note that this function considers 0 to be the best score. It will be inverted
// later.
void NaiveBayes::predict_column(pcolumn& ret, uint j) const {
    // Get the k-nearest columns
    proximity_queue q = d->knearest(j, k, (double)(max_distance));

    while (q.size() > 0) {
        id_dist_iter kth_j2 = q.top(); // get first  of the k-nearest columns

        // Get the observations in those columns
        size_t masked_j2_obs_size = kth_j2.matrix_iter->observations( kth_j2.id ).size();
        double masked_intersection_over_total = kth_j2.matrix_iter->intersection( j, kth_j2.id ).size() / double(masked_j2_obs_size);

        // Now that we've calculated the necessary values for the hypergeometric,
        // let's figure out the matrix we actually want to predict from, and get
        // those observations.
        id_set kth_j2_obs = predict_from( kth_j2.matrix_iter )->observations( kth_j2.id );

        // Iterate through the ACTUAL (unmasked) observations and predict the
        // score in applicable predict matrix cells using those observations and
        // the masked distance calculation.
        for (id_set::const_iterator it = kth_j2_obs.begin(); it != kth_j2_obs.end(); ++it) {
            pcolumn::iterator ret_it = ret.find(*it);

            // THIS IS THE ACTUAL NAIVE BAYES FORMULA!
            // (inside the Mult operator)
            float score_mod = (1.0 - masked_intersection_over_total * (1.0 - kth_j2.distance));

            if (ret_it == ret.end()) ret[*it] = score_mod;         // insert
            else                     ret_it->second *= score_mod;  // multiply
        }
        q.pop(); // move on to the next nearest item
    }
}

