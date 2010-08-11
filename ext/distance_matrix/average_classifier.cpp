#include "classifier.h"


// This is the NaiveBayes constructor as well.
AverageClassifier::AverageClassifier(const DistanceMatrix* const rhs, size_t k_, float max_distance_, float distance_exponent_)
        : Classifier(rhs), k(k_), max_distance(max_distance_), distance_exponent(distance_exponent_) { }


// Note that this function considers 0 to be the best score. It will be inverted
// later.
// Classifier does not take into account min_idf in its integration. This value
// only affects the actual "distance" calculation, not predictions.
void AverageClassifier::predict_column(pcolumn& ret, uint j) const {
    // Get the k-nearest columns
    proximity_queue q = d->knearest(j, k, (double)(max_distance));

    float q_size = (float)(q.size());

    while (q.size() > 0) {
        id_dist_iter kth_j2 = q.top(); // get first  of the k-nearest columns

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
            float score_mod = std::pow(kth_j2.distance, (double)(distance_exponent) );

            if (ret_it == ret.end()) ret[*it] = score_mod;         // insert
            else                     ret_it->second += score_mod;  // average
        }
        q.pop(); // move on to the next nearest item
    }

    // Average the scores
    for (pcolumn::iterator it = ret.begin(); it != ret.end(); ++it)
        it->second /= q_size;
}