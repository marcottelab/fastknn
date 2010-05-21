#include <map>
#include <set>
#include <algorithm>
using std::set;
using std::set_intersection;

#include "hypergeometric.h"
#include "euclidean.h"
#include "phenomatrix.h"
// typedef boost::numeric::ublas::mapped_matrix<double> dmatrix_t;

size_t tree_matrix_row_count(conn_t& c, uint id) {
    work_t w(c);
    return 0;
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


class DistanceMatrix {
public:
    DistanceMatrix(const string& dbstr, uint predict_matrix_id, uint source_matrix_id, string distfn = "hypergeometric")
    : c(new conn_t(dbstr)),
            // m(),
            source_matrix_(c, source_matrix_id, true),
            predict_matrix_(c, predict_matrix_id, false),
            distance_function(switch_distance_function(distfn))
    {
    }

    ~DistanceMatrix() {
        delete c;
    }

    //double distance(size_t k, size_t m, size_t n, size_t N) const {
    //    return (*distance_function)(k,m,n,N);
    //}

    // Return the distance, according to our distance function, between the two
    // columns.
    double distance(uint j1, uint j2) const {
        return (*distance_function)(
                predict_matrix_.observations_size(j1),
                source_matrix_.observations_size(j2),
                intersection_size(j1,j2),
                max_intersection_size());
    }

    // Find the single closest column in the source matrix to j in the predict matrix.
    std::pair<uint,double> nearest(uint j) const {
        id_set s = source_matrix_.column_ids();

        double min_dist = 100;
        uint min_dist_id = 0;

        for (id_set::const_iterator k = s.begin(); k != s.end(); ++k) {
            if (j == *k) continue; // Don't count it when the columns are the same
            double d_jk = distance(j, *k);
            if (d_jk < min_dist) {
                min_dist = d_jk;
                min_dist_id = *k;
            }
        }
        
        return std::make_pair<uint,double>(min_dist_id, min_dist);
    }

    // Find the k nearest columns
//    id_set knearest(uint j, size_t k) const {
//        id_set pj = observations(j);
//    }

    // Get the items that are common between j1 in predict matrix and j2 in
    // source matrix
    id_set intersection(uint j1, uint j2) const {
        id_set s1 = predict_matrix_.observations(j1);
        id_set s2 = source_matrix_.observations(j2);

        id_set ret;
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                         std::insert_iterator<id_set>(ret,ret.begin()));
        return ret;
    }

    // Count the number of items in common between j1 and j2
    size_t intersection_size(uint j1, uint j2) const {
        return intersection(j1, j2).size();
    }


#ifdef RICE
    uint rb_nearest(uint j) const {
        return nearest(j).first;
    }
    double rb_nearest_distance(uint j) const {
        return nearest(j).second;
    }
//    // These are the same as predict_parent_id and source_parent_id, but they
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

    size_t max_intersection_size() const { return source_matrix_.tree_row_count(); }
protected:
    conn_t* c;

    /// dmatrix_t m;

    // Database-contents-related stuff
    // These are protected because we don't want anyone accessing the connection
    // after this class has been destroyed. They'll have to create new ones.
    Phenomatrix source_matrix_;
    Phenomatrix predict_matrix_;

    // Allow different distance functions to be subbed in.
    double (*distance_function)(size_t, size_t, size_t, size_t);
};
