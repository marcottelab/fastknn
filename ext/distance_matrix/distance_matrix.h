#ifndef DISTANCE_MATRIX_H_
# define DISTANCE_MATRIX_H_


#include <map>
#include <set>
#include <algorithm>
using std::set;
using std::set_intersection;

typedef std::set<uint> id_set;

#include "phenomatrix_pair.h"
#include "cparams.h"
#include "classifier.h"
// typedef boost::numeric::ublas::mapped_matrix<double> dmatrix_t;


//size_t tree_matrix_row_count(conn_t& c, uint id) {
//    work_t w(c);
//    return 0;
//}


// Arguments for classifier:
// - k for knn
//float (*switch_classifier_function(const std::string& classifier))(size_t) {
//    std::map<std::string, float(*)(size_t)> choices;
//    choices["naivebayes"]     = &naivebayes;
//
//    return choices[classifier];
//}


class Classifier;

class DistanceMatrix {
    friend class Classifier;
    friend class NaiveBayes;
public:

    // This constructor allows a connection to be shared among multiple objects by
    // taking a pointer to an existing one. It is assumed that Ruby won't pass
    // such a pointer, and so we use a set of uints instead of a Rice::Array (which
    // is neccessary for the Ruby interface) as well.
    //
    // In other words, this constructor is exclusively for calling from within
    // a C++ environment of some sort.
    DistanceMatrix(uint, id_set, string, cparams);

    ~DistanceMatrix();

    // Find source matrix by columns
    matrix_list::const_iterator find_by_column(uint j) const {
        for (matrix_list::const_iterator dt = source_matrices.begin(); dt != source_matrices.end(); ++dt) {
#ifdef DEBUG_TRACE_INTERSECTION
            cerr << "distance_matrix.h: find_by_column: Checking matrix " << dt->id() << endl;
#endif
            if (dt->source_matrix_has_column(j)) {
#ifdef DEBUG_TRACE_INTERSECTION
                cerr << "\t...found on " << dt->id() << endl;
#endif
                return dt;
            }
        }
        cerr << "distance_matrix.h: Warning: source matrix with column " << j << " was not found." << endl;
        return source_matrices.end(); // not found
    }

    // Find source matrix by matrix ID
    matrix_list::const_iterator find(uint id) const {
        for (matrix_list::const_iterator dt = source_matrices.begin(); dt != source_matrices.end(); ++dt)
            if (dt->id() == id) return dt;
        cerr << "distance_matrix.h: Warning: source matrix with id " << id << " was not found." << endl;
        return source_matrices.end(); // not found
    }


    // Given a row i and a column j, calculate a score (using the classifier function)
    // for the gene's likelihood of being involved in the phenotype.
    pcolumn predict(uint j) const;


    // Return the distance, according to our distance function, between the two
    // columns.
    double distance(uint j1, uint j2) const {
        matrix_list::const_iterator j2source = find_by_column(j2);
        if (j2source == source_matrices.end()) {
            string err = "distance_matrix.h: distance(2): column supplied in second argument, " + lexical_cast<string>(j2) + ", does not appear to exist in any source matrix.";
#ifdef RICE
            throw Rice::Exception(rb_eArgError, err.c_str());
#else
            cerr << err << endl;
            throw;
#endif
        }

        return j2source->distance(j1, j2);
    }


    // Note that this function returns only the SINGLE NEAREST -- as in, the first
    // found!
    //
    // If you want multiple nearest, use a different function.
    id_dist_iter nearest(uint j) const {
        // Do a priming read to save an assignment -- particularly since we're
        // likely to have only one source matrix.
        matrix_list::const_iterator pt = source_matrices.begin();
        id_dist_iter min = pt->nearest(j);
        ++pt;

        for (; pt != source_matrices.end(); ++pt) {
            id_dist_iter min_tmp = pt->nearest(j);
            if (min_tmp.distance < min.distance) {
                min = min_tmp;
                min.matrix_iter = pt; // this keeps us from having to look up each matrix again
            }
        }

        return min;
    }

    // Get the items that are common between j1 in predict matrix and j2 in
    // source matrix
    id_set intersection(uint j1, uint j2) const {
        matrix_list::const_iterator f = find_by_column(j2);
#ifdef DEBUG_TRACE_INTERSECTION
        cerr << "intersection(2): source matrix = " << f->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
#endif
        if (f == source_matrices.end())
            return id_set(); // empty
        
        return f->intersection(j1, j2);
    }

    // Find the k items closest to j
    proximity_queue knearest(uint j, size_t k = 1, double kth_so_far = 100.0) const {
        proximity_queue q;

        for (matrix_list::const_iterator source_matrix_iter = source_matrices.begin(); source_matrix_iter != source_matrices.end(); ++source_matrix_iter) {
            // Add items on to the queue (q) that are within k (or kth_so_far)
            source_matrix_iter->knearest(q, j, k, kth_so_far, source_matrix_iter);
        }

        // This is the actual return-queue:
        proximity_queue ret;

        // Find the first k items
        while (k > 0) {
            kth_so_far = q.top().distance;
            ret.push(q.top());
            q.pop();
            k--;
        }

        // Now include further items equal to kth_so_far
        while (q.top().distance == kth_so_far) { ret.push(q.top()); q.pop(); }

        return ret;
    }
    
#ifdef RICE
    Rice::Object intersection_size(uint i, uint j) const {
        return to_ruby<size_t>(intersection(i, j).size());
    }

    Rice::Object source_matrix_ids() const {
        Rice::Array ids;
        for (matrix_list::const_iterator st = source_matrices.begin(); st != source_matrices.end(); ++st)
            ids.push(to_ruby<uint>(st->id()));
        return ids;
    }

#endif
protected:

    // Set up the classifier to use for predictions
    void construct_classifier(const cparams&);

    // Database-contents-related stuff
    matrix_list source_matrices;
    Phenomatrix predict_matrix_;

    // Allow different classifiers to be subbed in
    Classifier* classifier;
};


#endif // DISTANCE_MATRIX_H_
