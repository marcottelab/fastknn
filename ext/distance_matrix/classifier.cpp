#include "classifier.h"

Classifier::Classifier(const DistanceMatrix* const rhs, id_set row_ids)
: d(rhs)
{
    if (row_ids_.size() == 0)
        row_ids_ = d->predict_matrix_.row_ids();
}

boost::unordered_map<uint, float> Classifier::operator()(uint j) const {
    pcolumn ret;

    // Get the nearest column
    id_dist_iter nearest_j = d->nearest(j);

    // Get the observations in that column
    id_set nearest_j_obs = nearest_j.matrix_iter->observations( nearest_j.id );

    // Make a prediction
    for (id_set::const_iterator it = nearest_j_obs.begin(); it != nearest_j_obs.end(); ++it)
        ret[*it] = 0.5;

    // Fill in the rest with 0s.
    for (id_set::const_iterator it = row_ids_.begin(); it != row_ids_.end(); ++it)
        if (ret.find(*it) == ret.end()) ret[*it] = 0.0;

    return ret;
}


NaiveBayes::NaiveBayes(const DistanceMatrix* const rhs, uint k_, id_set row_ids)
: Classifier(rhs, row_ids), k(k_)
{ }

boost::unordered_map<uint,float> NaiveBayes::operator()(uint j) const {
    pcolumn ret;

    // Get the k-nearest columns
    proximity_queue q = d->knearest(j, k);

    while (q.size() > 0) {
        id_dist_iter j2 = q.top();

        // Get the observations in those columns
        id_set j2_obs = j2.matrix_iter->observations( j2.id );
        double intersection_over_total = d->intersection_size_given_matrix( j, j2.matrix_iter, j2.id ) / double(j2_obs.size());

        for (id_set::const_iterator it = j2_obs.begin(); it != j2_obs.end(); ++it) {
            pcolumn::iterator ret_it = ret.find(*it);

            // THIS IS THE ACTUAL NAIVE BAYES FORMULA!
            // (inside the Mult operator)
            float score_mod = (1.0 - intersection_over_total * (1.0 - j2.distance));

            if (ret_it == ret.end()) ret[*it] = score_mod;         // insert
            else                     ret_it->second *= score_mod;  // multiply
        }
        q.pop(); // move on to the next nearest item
    }

    return ret;
}