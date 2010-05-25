#ifndef CLASSIFIER_H_
# define CLASSIFIER_H_

#include "typedefs.h"
#include "distance_matrix.h"

class DistanceMatrix;

class Classifier {
public:
    // Can't change the DistanceMatrix or the passed-in pointer
    //
    // If an empty set is given, just predict all rows. Otherwise, only predict
    // the given rows.
    Classifier(const DistanceMatrix* const rhs, id_set row_ids = id_set());
    virtual ~Classifier() { }

    virtual boost::unordered_map<uint, float> operator()(uint j) const;
protected:

    const DistanceMatrix* d; // Can change the pointer but not the DistanceMatrix itself
    id_set row_ids_;
};

class NaiveBayes : public Classifier {
public:
    NaiveBayes(const DistanceMatrix* const rhs, uint k_, id_set row_ids = id_set());
    virtual ~NaiveBayes() { }

    virtual boost::unordered_map<uint, float> operator()(uint j) const;
protected:
    uint k;
};

#endif
