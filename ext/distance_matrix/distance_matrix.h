#ifndef DISTANCE_MATRIX_H_
# define DISTANCE_MATRIX_H_


#include <map>
#include <set>
#include <algorithm>
using std::set;
using std::set_intersection;

#include "typedefs.h"
#include "cparams.h"
#include "hypergeometric.h"
#include "euclidean.h"
#include "classifier.h"
// typedef boost::numeric::ublas::mapped_matrix<double> dmatrix_t;

typedef std::priority_queue<id_dist_iter, std::vector<id_dist_iter>, std::greater<vector<id_dist_iter>::value_type> > proximity_queue;


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
#ifdef RICE
    // This constructor is the one we use for the Ruby interface (RICE) since
    // Ruby would likely have trouble with std::set. Instead, it takes an array.
    DistanceMatrix(const string& dbstr, uint predict_matrix_id, const Rice::Array&, const string& distfn, Rice::Object);
#endif

    // This constructor allows a connection to be shared among multiple objects by
    // taking a pointer to an existing one. It is assumed that Ruby won't pass
    // such a pointer, and so we use a set of uints instead of a Rice::Array (which
    // is neccessary for the Ruby interface) as well.
    //
    // In other words, this constructor is exclusively for calling from within
    // a C++ environment of some sort.
    DistanceMatrix(conn_t* c_, uint predict_matrix_id, const id_set& source_matrix_ids, const string& distfn, const cparams& classifier_params);

    ~DistanceMatrix();

    //double distance(size_t k, size_t m, size_t n, size_t N) const {
    //    return (*distance_function)(k,m,n,N);
    //}

    // Find source matrix by columns
    matrix_list::const_iterator find_source_matrix_by_column(uint j) const {
        for (matrix_list::const_iterator dt = source_matrices.begin(); dt != source_matrices.end(); ++dt) {
#ifdef DEBUG_TRACE_INTERSECTION
            cerr << "distance_matrix.h: find_source_matrix_by_column: Checking matrix " << dt->id() << endl;
#endif
            if (dt->has_column(j)) {
#ifdef DEBUG_TRACE_INTERSECTION
                cerr << "\t...found on " << dt->id() << endl;
#endif
                return dt;
            }
        }
        cerr << "distance_matrix.h: Warning: source matrix with column " << j << " was not found." << endl;
        return source_matrices.end(); // not found
    }

    // Find source matrix by columns
    matrix_list::const_iterator find_source_matrix_by_id(uint id) const {
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
        matrix_list::const_iterator j2source = find_source_matrix_by_column(j2);
        if (j2source != source_matrices.end())
            return distance_given_matrix(j1, j2source, j2);
        return 1.0;
    }


    // Note that this function returns only the SINGLE NEAREST -- as in, the first
    // found!
    //
    // If you want multiple nearest, use a different function.
    id_dist_iter nearest(uint j) const {
        // Do a priming read to save an assignment -- particularly since we're
        // likely to have only one source matrix.
        matrix_list::const_iterator pt = source_matrices.begin();
        id_dist_iter min = nearest_given_matrix(pt, j);
        ++pt;

        for (; pt != source_matrices.end(); ++pt) {
            id_dist_iter min_tmp = nearest_given_matrix(pt, j);
            if (min_tmp.distance < min.distance)
                min = min_tmp;
        }

        return min;
    }

    // Get the items that are common between j1 in predict matrix and j2 in
    // source matrix
    id_set intersection(uint j1, uint j2) const {
        matrix_list::const_iterator f = find_source_matrix_by_column(j2);
#ifdef DEBUG_TRACE_INTERSECTION
        cerr << "intersection(2): source matrix = " << f->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
#endif
        if (f == source_matrices.end())
            return id_set(); // empty
        
        return intersection_given_matrix(j1, f, j2);
    }

    // Count the number of items in common between j1 and j2 (j1 in predict matrix, j2 in a source matrix)
    size_t intersection_size(uint j1, uint j2) const {
        matrix_list::const_iterator f = find_source_matrix_by_column(j2);
#ifdef DEBUG_TRACE_INTERSECTION
        cerr << "intersection_size(2): source matrix = " << f->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
#endif
        
        if (f == source_matrices.end())
            return 0; // empty

        return intersection_size_given_matrix(j1, f, j2);
    }


#ifdef RICE
    Rice::Object rb_predict(uint j) const {
        pcolumn predictions = predict(j);
        Rice::Hash h;
        for (pcolumn::const_iterator i = predictions.begin(); i != predictions.end(); ++i)
            h[ to_ruby<uint>(i->first) ] = to_ruby<float>(i->second);
        return h;
    }


    Rice::Object rb_knearest(uint j, size_t k = 1) const {
        proximity_queue q = knearest(j, k);
        Rice::Array ret;
        while (q.size() > 0) {
            id_dist_iter di = q.top();
            ret.push( di.to_a() );
            q.pop();
        }

        return ret;
    }
    
    // Returns an id or nil
    Rice::Object rb_nearest_id(uint j) const {
        id_dist_iter n = nearest(j);
        return (n.id == 0) ? Rice::Object() : Rice::Object( to_ruby<uint>(n.id) );
    }
    // Returns a double or nil
    Rice::Object rb_nearest_distance(uint j) const {
        id_dist_iter n = nearest(j);
        return (n.id == 0) ? Rice::Object() : Rice::Object( to_ruby<double>(n.distance) );
    }

    // Returns both the id and the distance, or nil
    Rice::Object rb_nearest(uint j) const {
        id_dist_iter n = nearest(j);
        return (n.id == 0) ? Rice::Object() : n.to_a();
    }

    // Return the distance, according to our distance function, between the two
    // columns.
    Rice::Object rb_distance(uint j1, uint j2) const {
        matrix_list::const_iterator j2source = find_source_matrix_by_column(j2);
        if (j2source != source_matrices.end())
            return to_ruby<double>((*distance_function)(
                                        predict_matrix_.observations_size(j1),
                                        j2source->observations_size(j2),
                                        intersection_size_given_matrix(j1, j2source, j2),
                                        max_intersection_size_given_matrix(j2source)));
        return Rice::Object(); // nil
    }

    Rice::Array source_matrix_ids() const {
        Rice::Array ids;
        for (matrix_list::const_iterator st = source_matrices.begin(); st != source_matrices.end(); ++st)
            ids.push(to_ruby<uint>(st->id()));
        return ids;
    }

//    // These are the same as predict_parent_id and source_parent_id, but theyfind_
//    // return Ruby nil instead of 0
//    Rice::Object rb_predict_parent_id() {
//        uint id = predict_parent_id();
//        return (id == 0) ? Rice::Object() : Rice::Object(to_ruby<uint>(id));
//    }
//    Rice::Object rb_source_parent_id() {
//        uint id = source_parent_id();
//        return (id == 0) ? Rice::Object() : Rice::Object(to_ruby<uint>(id));
//    }
//    Rice::Object rb_predict_root_id() {
//        uint id = predict_root_id();
//        return (id == 0) ? Rice::Object() : Rice::Object(to_ruby<uint>(id));
//    }
//    Rice::Object rb_source_root_id() {
//        uint id = source_root_id();
//        return (id == 0) ? Rice::Object() : Rice::Object(to_ruby<uint>(id));
//    }
#endif

    //Phenomatrix source_matrix() const { return source_matrix_; }
    //Phenomatrix predict_matrix() const { return predict_matrix_; }


    size_t max_intersection_size() const {
        size_t max_size = 0;
        
        for (matrix_list::const_iterator pt = source_matrices.begin(); pt != source_matrices.end(); ++pt) {
            size_t tmp = max_intersection_size_given_matrix(pt);
            if (tmp > max_size)
                max_size = tmp;
        }
        return max_size;
    }
    
protected:

    // Set up the classifier to use for predictions
    void construct_classifier(const cparams&);


    double distance_given_matrix(const uint& j1, matrix_list::const_iterator source_matrix_iter, const uint& j2) const {
        return (*distance_function)(
                    predict_matrix_.observations_size(j1),
                    source_matrix_iter->observations_size(j2),
                    intersection_size_given_matrix(j1,source_matrix_iter, j2),
                    max_intersection_size_given_matrix(source_matrix_iter));
    }


    // Find the single closest column in the source matrix to j in the predict matrix.
    id_dist_iter nearest_given_matrix(matrix_list::const_iterator source_matrix_iter, const uint& j) const {
        id_set s = source_matrix_iter->column_ids();

        double min_dist = 100;
        uint min_dist_id = 0;

        for (id_set::const_iterator k = s.begin(); k != s.end(); ++k) {
            if (j == *k) continue; // Don't count it when the columns are the same
            double d_jk = distance_given_matrix(j, source_matrix_iter, *k);
            if (d_jk < min_dist) {
                min_dist = d_jk;
                min_dist_id = *k;
            }
        }

        return id_dist_iter(min_dist_id, min_dist, source_matrix_iter);
    }


    // Find the k items closest to j
    proximity_queue knearest(const uint& j, size_t k = 1, double kth_so_far = 100.0) const {
        proximity_queue q;
        
        for (matrix_list::const_iterator source_matrix_iter = source_matrices.begin(); source_matrix_iter != source_matrices.end(); ++source_matrix_iter) {
            // Add items on to the queue (q) that are within k (or kth_so_far)
            knearest_given_matrix(q, source_matrix_iter, j, k, kth_so_far);
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


    // Find the k closest columns, and those of equivalent distance, to j
    // In other words, if we get to the kth closest and there's another column
    // with the same distance, we include both.
    void knearest_given_matrix(proximity_queue& q, matrix_list::const_iterator source_matrix_iter, const uint& j, const size_t& k, double& kth_so_far) const {
        id_set s = source_matrix_iter->column_ids();

        uint current_k = 0;

        // Only sort all on the first matrix. After that, can limit to those within
        // kth_so_far.
        for (id_set::const_iterator st = s.begin(); st != s.end(); ++st) {
            if (j == *st) continue; // Don't count it when the columns are the same
            double d_jk = distance_given_matrix(j, source_matrix_iter, *st);

            if (d_jk <= kth_so_far) { // don't add distances that are already outside the range
                q.push( id_dist_iter(*st, d_jk, source_matrix_iter) );

                // We can set a boundary within k elements, and not accept anymore
                // outside of it.
                ++current_k;
                if (current_k == k) kth_so_far = d_jk;
            }
        }
    }


    size_t max_intersection_size_given_matrix(matrix_list::const_iterator source_matrix_iter) const {
        return source_matrix_iter->tree_row_count();
    }
    

    id_set intersection_given_matrix(const uint& j1, matrix_list::const_iterator source_matrix_iter, const uint& j2) const {
#ifdef DEBUG_TRACE_INTERSECTION
        cerr << "distance_matrix.h: intersection_given_matrix(3): source matrix = " << source_matrix_iter->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
#endif

        id_set s1 = predict_matrix_.observations(j1);
        id_set s2 = source_matrix_iter->observations(j2);

        id_set ret;
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                         std::insert_iterator<id_set>(ret,ret.begin()));
        return ret;
    }

    size_t intersection_size_given_matrix(const uint& j1, matrix_list::const_iterator source_matrix_iter, const uint& j2) const {
#ifdef DEBUG_TRACE_INTERSECTION
        cerr << "distance_matrix.h: intersection_size_given_matrix(3): source matrix = " << source_matrix_iter->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
#endif
        return intersection_given_matrix(j1, source_matrix_iter, j2).size();
    }

    // Returns a function pointer to a distance function based on a request made via
    // a string.
    double (*switch_distance_function(const std::string& distance_measure))(size_t,size_t,size_t,size_t) {
        std::map<std::string, double(*)(size_t,size_t,size_t,size_t)> choices;
        choices["hypergeometric"] = &hypergeometric;
        choices["euclidean"]      = &euclidean;
        choices["manhattan"]      = &manhattan;

        return choices[distance_measure];
    }


    conn_t* c;
    bool destroy_c; // keep the connection intact upon destruction?
    /// dmatrix_t m;

    // Database-contents-related stuff
    // These are protected because we don't want anyone accessing the connection
    // after this class has been destroyed. They'll have to create new ones.
    matrix_list source_matrices;
    Phenomatrix predict_matrix_;

    // Allow different distance functions to be subbed in.
    double (*distance_function)(size_t, size_t, size_t, size_t);

    // Allow different classifiers to be subbed in
    Classifier* classifier;
};


#endif // DISTANCE_MATRIX_H_
