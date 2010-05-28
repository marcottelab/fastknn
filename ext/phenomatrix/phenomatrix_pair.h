#ifndef PHENOMATRIX_PAIR_H_
# define PHENOMATRIX_PAIR_H_

#include <stack>
using std::stack;

#include "id_dist.h"
#include "connection.h"
#include "typedefs.h"
#include "hypergeometric.h"
#include "euclidean.h"
#include "phenomatrix.h"


class id_dist;
class id_dist_iter;
typedef std::priority_queue<id_dist_iter, std::vector<id_dist_iter>, std::greater<vector<id_dist_iter>::value_type> > proximity_queue;

#ifndef MATRIX_LIST_DEFINED
# define MATRIX_LIST_DEFINED
class PhenomatrixPair;
typedef std::list<PhenomatrixPair>          matrix_list;
#endif

// Handles creation of phenomatrices for pairs of species such that they have the
// correct number of rows.
//
// Also allows "temporary" removal of rows for the purpose of calculating
// distances and such (see push_mask and pop_mask).
class PhenomatrixPair {
public:
    PhenomatrixPair(uint id, uint given_id, const string& distance_fn)
    : p(create_phenomatrix(id, given_id)),      // predict species matrix
      s(create_phenomatrix(given_id, given_id)),// source species matrix
      distance_function(switch_distance_function(distance_fn))
    { }

    PhenomatrixPair(const string& dbstr, uint id, uint given_id, const string& distance_fn)
    : p(create_phenomatrix(id, given_id)),
      s(create_phenomatrix(given_id, given_id)),
      distance_function(switch_distance_function(distance_fn))
    { }
    
    // Return the ID of the source matrix
    uint id() const {
        return s.top().id();
    }


    // Removes rows from the matrices on which we're calculating distances.
    void push_mask(id_set mask_rows) {
        p.push(Phenomatrix(p.top(), mask_rows));
        s.push(Phenomatrix(s.top(), mask_rows));
    }


    // Restores removed rows from the matrices on which we're calculating distances.
    bool pop_mask() {
        if (p.size() == 1) return false;
        p.pop();
        s.pop();
        return true;
    }


    // Return the distance, according to our distance function, between the two
    // columns.
    double distance(uint j1, uint j2) const {
        size_t obs_j1 = p.top().observations_size(j1);
        size_t obs_j2 = s.top().observations_size(j2);

        if (obs_j1 == 0 || obs_j2 == 0) return MAX_DISTANCE;

        return (*distance_function)( obs_j1, obs_j2, intersection(j1, j2).size(), max_intersection_size() );
    }

    bool source_matrix_has_column(uint j) const {
        return s.top().has_column(j);
    }

    bool predict_matrix_has_column(uint j) const {
        return p.top().has_column(j);
    }

    // Find the single nearest neighbor (FIRST FOUND, not equivalence class)
    id_dist_iter nearest(uint j) const;


    // This is more of a helper for knearest. It needs to be set up properly in
    // DistanceMatrix.
    void knearest(proximity_queue& q, const uint& j, const size_t& k, double& kth_so_far, matrix_list::const_iterator this_iter) const;


    id_set intersection(uint j1, uint j2) const {
        const id_set& s1 = p.top().observations(j1);
        const id_set& s2 = s.top().observations(j2);

        id_set ret;
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                          std::insert_iterator<id_set>(ret,ret.begin()));
        return ret;
    }

    size_t max_intersection_size() const {
        return s.top().row_count();
    }

    id_set observations(uint j) const {
        return s.top().observations(j);
    }
protected:

    // Returns a function pointer to a distance function based on a request made via
    // a string.
    double (*switch_distance_function(const std::string& distance_measure))(size_t,size_t,size_t,size_t) {
        std::map<std::string, double(*)(size_t,size_t,size_t,size_t)> choices;
        choices["hypergeometric"] = &hypergeometric;
        choices["euclidean"]      = &euclidean;
        choices["manhattan"]      = &manhattan;

        return choices[distance_measure];
    }

    static stack<Phenomatrix> create_phenomatrix(uint id, uint given_id) {
        stack<Phenomatrix> new_stack; new_stack.push(Phenomatrix(id, given_id));
        return new_stack;
    }

    stack<Phenomatrix> p;
    stack<Phenomatrix> s;

    // Allow different distance functions to be subbed in.
    double (*distance_function)(size_t, size_t, size_t, size_t);
};


#ifndef MATRIX_LIST_DEFINED
# define MATRIX_LIST_DEFINED
typedef std::list<PhenomatrixPair>          matrix_list;
#endif

#endif