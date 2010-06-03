#include "classifier.h"

void SimpleClassifier::predict_column(pcolumn& ret, uint j) const {
    // Get the nearest column
    id_dist_iter nearest_j = d->nearest(j);

    // Get the observations in that column
    id_set nearest_j_obs = predict_from(nearest_j.matrix_iter)->observations( nearest_j.id );

    // Make a prediction
    for (id_set::const_iterator it = nearest_j_obs.begin(); it != nearest_j_obs.end(); ++it)
        ret[*it] = 0.5;
}