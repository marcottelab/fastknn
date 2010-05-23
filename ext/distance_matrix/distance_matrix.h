//#define RICE

#ifdef RICE
#include <rice/Object.hpp>
#include <rice/Array.hpp>
using Rice::Array;
using Rice::Object;
#endif

#include <list>
#include <map>
#include <set>
#include <algorithm>
using std::set;
using std::list;
using std::set_intersection;

#include "hypergeometric.h"
#include "euclidean.h"
#include "../phenomatrix/phenomatrix.h"
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
#ifdef RICE
    DistanceMatrix(const string& dbstr, uint predict_matrix_id, const Array& source_matrix_ids, const string& distfn)
    : c(new conn_t(dbstr)), destroy_c(true), predict_matrix_(c, predict_matrix_id, false),
            distance_function(switch_distance_function(distfn))
    {
        for (Array::const_iterator st = source_matrix_ids.begin(); st != source_matrix_ids.end(); ++st) {
            uint id = from_ruby<uint>(*st);
            cerr << "distance_matrix.h: Adding phenomatrix " << id << " to distance matrix" << endl;
            source_matrices.push_back( Phenomatrix(c, id, true) );
        }
    }
#else
    DistanceMatrix(conn_t* c_, uint predict_matrix_id, const id_set& source_matrix_ids, const string& distfn = "hypergeometric")
            : c(c_), destroy_c(false), predict_matrix_(c, predict_matrix_id, false),
            distance_function(switch_distance_function(distfn))
    {
        for (id_set::const_iterator st = source_matrix_ids.begin(); st != source_matrix_ids.end(); ++st) {
            source_matrices.push_back( Phenomatrix(c, *st, true) );
        }
    }
#endif

    ~DistanceMatrix() {
        if (destroy_c) delete c;
    }

    //double distance(size_t k, size_t m, size_t n, size_t N) const {
    //    return (*distance_function)(k,m,n,N);
    //}

    // Find source matrix by columns
    list<Phenomatrix>::const_iterator find_source_matrix_by_column(uint j) const {
        for (list<Phenomatrix>::const_iterator dt = source_matrices.begin(); dt != source_matrices.end(); ++dt) {
            cerr << "distance_matrix.h: find_source_matrix_by_column: Checking matrix " << dt->id() << endl;
            if (dt->has_column(j)) {
                cerr << "\t...found on " << dt->id() << endl;
                return dt;
            }
        }
        cerr << "distance_matrix.h: Warning: source matrix with column " << j << " was not found." << endl;
        return source_matrices.end(); // not found
    }

    // Find source matrix by columns
    list<Phenomatrix>::const_iterator find_source_matrix_by_id(uint id) const {
        for (list<Phenomatrix>::const_iterator dt = source_matrices.begin(); dt != source_matrices.end(); ++dt)
            if (dt->id() == id) return dt;
        cerr << "distance_matrix.h: Warning: source matrix with id " << id << " was not found." << endl;
        return source_matrices.end(); // not found
    }

    // Return the distance, according to our distance function, between the two
    // columns.
    double distance(uint j1, uint j2) const {
        list<Phenomatrix>::const_iterator j2source = find_source_matrix_by_column(j2);
        if (j2source != source_matrices.end())
            return distance_given_matrix(j1, j2source, j2);
        return 1.0;
    }

    double distance_given_matrix(const uint& j1, list<Phenomatrix>::const_iterator source_matrix_iter, const uint& j2) const {
        return (*distance_function)(
                    predict_matrix_.observations_size(j1),
                    source_matrix_iter->observations_size(j2),
                    intersection_size_given_matrix(j1,source_matrix_iter, j2),
                    max_intersection_size_given_matrix(source_matrix_iter));
    }


    // Find the single closest column in the source matrix to j in the predict matrix.
    std::pair<uint,double> nearest_given_matrix(list<Phenomatrix>::const_iterator source_matrix_iter, const uint& j) const {
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
        
        return std::make_pair<uint,double>(min_dist_id, min_dist);
    }


    std::pair<uint,double> nearest(uint j) const {
        std::pair<uint,double> min;

        for (list<Phenomatrix>::const_iterator pt = source_matrices.begin(); pt != source_matrices.end(); ++pt) {
            std::pair<uint,double> min_tmp = nearest_given_matrix(pt, j);
            if (min_tmp.second < min.second) min = min_tmp;
        }

        return min;
    }


    // Find the k nearest columns
//    id_set knearest(uint j, size_t k) const {
//        id_set pj = observations(j);
//    }


    // Get the items that are common between j1 in predict matrix and j2 in
    // source matrix
    id_set intersection(uint j1, uint j2) const {
        list<Phenomatrix>::const_iterator f = find_source_matrix_by_column(j2);
        cerr << "intersection(2): source matrix = " << f->id() << ", j1=" << j1 << ", j2=" << j2 << endl;

        if (f == source_matrices.end())
            return id_set(); // empty
        
        return intersection_given_matrix(j1, f, j2);
    }

    // Count the number of items in common between j1 and j2 (j1 in predict matrix, j2 in a source matrix)
    size_t intersection_size(uint j1, uint j2) const {
        list<Phenomatrix>::const_iterator f = find_source_matrix_by_column(j2);
        cerr << "intersection_size(2): source matrix = " << f->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
        
        if (f == source_matrices.end())
            return 0; // empty

        return intersection_size_given_matrix(j1, f, j2);
    }


#ifdef RICE
    uint rb_nearest(uint j) const {
        return nearest(j).first;
    }
    double rb_nearest_distance(uint j) const {
        return nearest(j).second;
    }

    // Return the distance, according to our distance function, between the two
    // columns.
    Rice::Object rb_distance(uint j1, uint j2) const {
        list<Phenomatrix>::const_iterator j2source = find_source_matrix_by_column(j2);
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
        for (list<Phenomatrix>::const_iterator st = source_matrices.begin(); st != source_matrices.end(); ++st)
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
        
        for (list<Phenomatrix>::const_iterator pt = source_matrices.begin(); pt != source_matrices.end(); ++pt) {
            size_t tmp = max_intersection_size_given_matrix(pt);
            if (tmp > max_size)
                max_size = tmp;
        }
        return max_size;
    }
    
protected:
    size_t max_intersection_size_given_matrix(list<Phenomatrix>::const_iterator source_matrix_iter) const {
        return source_matrix_iter->tree_row_count();
    }

    id_set intersection_given_matrix(const uint& j1, list<Phenomatrix>::const_iterator source_matrix_iter, const uint& j2) const {
        cerr << "distance_matrix.h: intersection_given_matrix(3): source matrix = " << source_matrix_iter->id() << ", j1=" << j1 << ", j2=" << j2 << endl;

        id_set s1 = predict_matrix_.observations(j1);
        id_set s2 = source_matrix_iter->observations(j2);

        id_set ret;
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                         std::insert_iterator<id_set>(ret,ret.begin()));
        return ret;
    }

    size_t intersection_size_given_matrix(const uint& j1, list<Phenomatrix>::const_iterator source_matrix_iter, const uint& j2) const {
        cerr << "distance_matrix.h: intersection_size_given_matrix(3): source matrix = " << source_matrix_iter->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
        return intersection_given_matrix(j1, source_matrix_iter, j2).size();
    }


    conn_t* c;
    bool destroy_c; // keep the connection intact upon destruction?
    /// dmatrix_t m;

    // Database-contents-related stuff
    // These are protected because we don't want anyone accessing the connection
    // after this class has been destroyed. They'll have to create new ones.
    list<Phenomatrix> source_matrices;
    Phenomatrix predict_matrix_;

    // Allow different distance functions to be subbed in.
    double (*distance_function)(size_t, size_t, size_t, size_t);
};
