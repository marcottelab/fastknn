#include "classifier.h"

Classifier::Classifier(const DistanceMatrix* const rhs)
: d(rhs)
{
    if (row_ids_.size() == 0)
        row_ids_ = d->predict_matrix_.row_ids();
}

