#include "classifier.h"


// This is the NaiveBayes constructor as well.
AverageClassifier::AverageClassifier(const DistanceMatrix* const rhs, size_t k_, float max_distance_)
        : Classifier(rhs), k(k_), max_distance(max_distance_) { }


// Note that this function considers 0 to be the best score. It will be inverted
// later.
// Classifier does not take into account min_idf in its integration. This value
// only affects the actual "distance" calculation, not predictions.
void AverageClassifier::predict_column(pcolumn& ret, uint j) const {
    typedef map<uint,float> opcolumn; // ordered phenotype column

    // Get the k-nearest columns
    proximity_queue q = d->knearest(j, k, (double)(max_distance));

    float q_size = (float)(q.size());

    opcolumn oret; // ordered map of return values


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
            opcolumn::iterator oret_it = oret.find(*it);

            if (oret_it == oret.end()) oret[*it] = 1;         // insert
            else                       oret_it->second += 1;  // average
        }
        q.pop(); // move on to the next nearest item
    }

    // Average the scores and invert them, then insert them in the return unordered_map
    for (opcolumn::iterator it = oret.begin(); it != oret.end(); ++it)
        ret[it->first] = 1.0 - (it->second / q_size);
}