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
    Classifier(const DistanceMatrix* const rhs);
    virtual ~Classifier() { }

    boost::unordered_map<uint, float> operator()(uint j) const {
        pcolumn ret;

        predict_column(ret, j);

        // Fill in the rest with 0s.
        zero_fill_remainder(ret);

        return ret;
    }
protected:
    // This is the function that needs to be overridden by other types of classifiers.
    virtual void predict_column(pcolumn& ret, uint j) const;

    list<Phenomatrix>::const_reverse_iterator predict_from(matrix_list::const_iterator source_pair) const {
        // If there is only one element, there is only one thing to predict.
        if (source_pair->size() == 1)
            return list<Phenomatrix>::const_reverse_iterator( source_pair->s_rbegin() );

        // If there are more elements, we don't want to predict the top -- we just
        // want to use that for distances and such.
        list<Phenomatrix>::const_reverse_iterator iter( source_pair->s_rbegin() );

        ++iter;
        return iter;
    }

    void zero_fill_remainder(pcolumn& ret) const {
        // Fill in the rest with 0s.
        for (id_set::const_iterator it = row_ids_.begin(); it != row_ids_.end(); ++it)
            if (ret.find(*it) == ret.end()) ret[*it] = 0.0;
    }


    const DistanceMatrix* d; // Can change the pointer but not the DistanceMatrix itself
    id_set row_ids_;
};


class NaiveBayes : public Classifier {
public:
    NaiveBayes(const DistanceMatrix* const rhs, size_t k_);

    virtual ~NaiveBayes() { }

protected:
    virtual void predict_column(pcolumn& ret, uint j) const;

    size_t k;
};


#endif
