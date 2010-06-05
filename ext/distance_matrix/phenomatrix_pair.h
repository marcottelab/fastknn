#ifndef PHENOMATRIX_PAIR_H_
# define PHENOMATRIX_PAIR_H_

#include <stack>
#include <list>
using std::list;
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
    : p(create_phenomatrix_stack(id, given_id)),      // predict species matrix
      s(create_phenomatrix_list(given_id, given_id)), // source species matrix
      distance_function(switch_distance_function(distance_fn))
    { }

    PhenomatrixPair(const string& dbstr, uint id, uint given_id, const string& distance_fn)
    : p(create_phenomatrix_stack(id, given_id)),
      s(create_phenomatrix_list(given_id, given_id)),
      distance_function(switch_distance_function(distance_fn))
    { }

    // Copy constructor
    PhenomatrixPair(const PhenomatrixPair& rhs)
    : p(rhs.p), s(rhs.s), distance_function(rhs.distance_function) { }
    
    // Return the ID of the source matrix
    uint id() const {
        return s.back().id();
    }


    // Removes rows from the matrices on which we're calculating distances.
    void push_mask(id_set mask_rows) {
        p.push(      Phenomatrix( p.top() , mask_rows ) );
        s.push_back( Phenomatrix( s.back(), mask_rows ) );
    }


    // Restores removed rows from the matrices on which we're calculating distances.
    bool pop_mask() {
        if (p.size() == 1) return false;
        p.pop();
        s.pop_back();
        return true;
    }


    // Return the distance, according to our distance function, between the two
    // columns.
    double distance(uint j1, uint j2) const {
        size_t obs_j1 = p.top().observations_size(j1);
        size_t obs_j2 = s.back().observations_size(j2);

        if (obs_j1 == 0 || obs_j2 == 0) return MAX_DISTANCE;

        return (*distance_function)( obs_j1, obs_j2, intersection(j1, j2).size(), max_intersection_size() );
    }

    bool source_matrix_has_column(uint j) const {
        return s.back().has_column(j);
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
        const id_set& s2 = s.back().observations(j2);

        id_set ret;
        set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                          std::insert_iterator<id_set>(ret,ret.begin()));
        return ret;
    }

    size_t max_intersection_size() const {
        return s.back().row_count();
    }

    id_set observations(uint j) const {
        return s.back().observations(j);
    }

    list<Phenomatrix>::const_reverse_iterator s_rbegin() const {
        return s.rbegin();
    }

    // Return the number of matrices on the stack.
    size_t size() const {
        return s.size();
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

    static stack<Phenomatrix> create_phenomatrix_stack(uint id, uint given_id) {
        stack<Phenomatrix> new_stack; new_stack.push(Phenomatrix(id, given_id));
        return new_stack;
    }

    static list<Phenomatrix> create_phenomatrix_list(uint id, uint given_id) {
        list<Phenomatrix> new_list; new_list.push_back(Phenomatrix(id, given_id));
        return new_list;
    }

    stack<Phenomatrix> p;
    list<Phenomatrix> s;

    // Allow different distance functions to be subbed in.
    double (*distance_function)(size_t, size_t, size_t, size_t);
};


#ifndef MATRIX_LIST_DEFINED
# define MATRIX_LIST_DEFINED
typedef std::list<PhenomatrixPair>          matrix_list;
#endif

#endif