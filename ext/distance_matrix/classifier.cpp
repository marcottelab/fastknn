#include "classifier.h"

Classifier::Classifier(const DistanceMatrix* const rhs)
: d(rhs)
{
    if (row_ids_.size() == 0)
        row_ids_ = d->predict_matrix_.row_ids();
}

void Classifier::predict_column(pcolumn& ret, uint j) const {
    // Get the nearest column
    id_dist_iter nearest_j = d->nearest(j);

    // Get the observations in that column
    id_set nearest_j_obs = predict_from(nearest_j.matrix_iter)->observations( nearest_j.id );

    // Make a prediction
    for (id_set::const_iterator it = nearest_j_obs.begin(); it != nearest_j_obs.end(); ++it)
        ret[*it] = 0.5;
}
