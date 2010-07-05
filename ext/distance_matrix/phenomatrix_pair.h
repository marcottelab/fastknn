#ifndef PHENOMATRIX_PAIR_H_
# define PHENOMATRIX_PAIR_H_

#ifdef RICE
#include <rice/Object.hpp>
#include <rice/Symbol.hpp>
#include <rice/Data_Type.hpp>
#endif

#include <iostream>
#include <stack>
#include <list>
#include <boost/foreach.hpp>
using std::list;
using std::stack;
using std::cerr;
using std::endl;
using std::pair;
using std::make_pair;

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
#ifdef RICE
    // Specialty constructor just for Ruby. Allows Ruby to pass in
    PhenomatrixPair(Rice::Object self, Rice::Object predict_or_id, Rice::Object source_or_id, size_t min_genes)
    : s(create_phenomatrix_list(source_or_id, min_genes)),
      p(create_phenomatrix_stack(predict_or_id, s.back().id(), min_genes)),
      distance_function(switch_distance_function("hypergeometric"))
    { }
#endif
    PhenomatrixPair(uint id, PhenomatrixBase given_phenomatrix, size_t min_genes = 2)
    : s(create_phenomatrix_list(given_phenomatrix)),
      p(create_phenomatrix_stack(id, given_phenomatrix.id(), min_genes)),
      distance_function(switch_distance_function("hypergeometric"))
    { }

    PhenomatrixPair(uint id, uint given_id, size_t min_genes = 2)
    : s(create_phenomatrix_list(given_id, min_genes)),            // source species matrix
      p(create_phenomatrix_stack(id, given_id, min_genes)),       // predict species matrix
      distance_function(switch_distance_function("hypergeometric"))
    { }

    // Copy constructor
    PhenomatrixPair(const PhenomatrixPair& rhs)
    : s(rhs.s),
      p(rhs.p),
      distance_function(rhs.distance_function)
    { }

    ~PhenomatrixPair() { }
    
    // Return the ID of the source matrix
    uint id() const {
        return s.back().id();
    }


    // Removes rows from the matrices on which we're calculating distances.
    void push_mask(id_set mask_rows) {
        p.push(      Phenomatrix(     p.top() , mask_rows ) );
        s.push_back( PhenomatrixBase( s.back(), mask_rows ) );
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

    pair<size_t, size_t> min_observations_count() const {
        return make_pair<size_t,size_t>(p.top().min_observations_count(), s.back().min_observations_count());
    }

    id_set observations(uint j) const {
        return s.back().observations(j);
    }

    list<PhenomatrixBase>::const_reverse_iterator s_rbegin() const {
        return s.rbegin();
    }

    // Return the number of matrices on the stack.
    size_t size() const {
        return s.size();
    }

    id_set row_ids() const {
        return s.back().row_ids();
    }

    id_set predictable_column_ids() const {
        return p.top().column_ids();
    }

#ifdef RICE
    Rice::Object get_distance_function() const {
        std::map<double(*)(size_t,size_t,size_t,size_t), std::string> choices;
        choices[&hypergeometric] = "hypergeometric";
        choices[&euclidean]      = "euclidean";
        choices[&manhattan]      = "manhattan";

        return Rice::Symbol(choices[distance_function]);
    }

    // Take a symbol, e.g., :hypergeometric, and use it to set the distance function.
    void set_distance_function(Rice::Object dfn) {
        distance_function = switch_distance_function( from_ruby<Rice::Symbol>(dfn).str() );
    }
#endif
    
    // Take a string, e.g., "hypergeometric", and use it to set the distance function.
    void set_distance_function_str(const string& dfn) {
        distance_function = switch_distance_function(dfn);
    }
protected:

    // Returns a function pointer to a distance function based on a request made via
    // a string.
    static double (*switch_distance_function(const std::string& distance_measure))(size_t,size_t,size_t,size_t) {
        std::map<std::string, double(*)(size_t,size_t,size_t,size_t)> choices;
        choices["hypergeometric"] = &hypergeometric;
        choices["euclidean"]      = &euclidean;
        choices["manhattan"]      = &manhattan;

        return choices[distance_measure];
    }

    // Create the predict matrix stack
    static stack<Phenomatrix> create_phenomatrix_stack(uint id, uint given_id, size_t min_genes) {
        stack<Phenomatrix> new_stack; new_stack.push(Phenomatrix(id, given_id));
        return new_stack;
    }

    // Create the source matrix stack (well, actually a list)
    static list<PhenomatrixBase> create_phenomatrix_list(uint id, size_t min_genes) {
        list<PhenomatrixBase> new_list; new_list.push_back(PhenomatrixBase(id));
        return new_list;
    }
    static list<PhenomatrixBase> create_phenomatrix_list(const PhenomatrixBase& pb) {
        list<PhenomatrixBase> new_list; new_list.push_back(pb);
        return new_list;
    }
    
#ifdef RICE
    static stack<Phenomatrix> create_phenomatrix_stack(Rice::Object predict_or_id, uint given_id, size_t min_genes) {
        if (predict_or_id.is_a( rb_cFixnum ))
            return create_phenomatrix_stack( from_ruby<uint>(predict_or_id), given_id, min_genes );

        else if (predict_or_id.is_a( Rice::Data_Type<Phenomatrix>::klass() )) {
            stack<Phenomatrix> new_stack;
            new_stack.push( from_ruby<Phenomatrix>(predict_or_id) );

            // Check that the predict matrix has the source matrix as a given
            if (new_stack.top().source_id() != given_id)
                throw Rice::Exception(rb_eArgError, "phenomatrix_pair.h: create_phenomatrix_stack: Source matrix does not match the source id given to the predict matrix");

            return new_stack;

        } else
          throw Rice::Exception(rb_eArgError, "phenomatrix_pair.h: create_phenomatrix_stack: Argument must be a Phenomatrix");  
    }

    static list<PhenomatrixBase> create_phenomatrix_list(Rice::Object source_or_id, size_t min_genes) {
        if (source_or_id.is_a( rb_cFixnum ))
            return create_phenomatrix_list(from_ruby<uint>(source_or_id), min_genes);

        else if (source_or_id.is_instance_of(Rice::Data_Type<PhenomatrixBase>::klass())) {
            list<PhenomatrixBase> new_list;
            new_list.push_back(from_ruby<PhenomatrixBase>(source_or_id));
            return new_list;
            
        } else
          throw Rice::Exception(rb_eArgError, "phenomatrix_pair.h: create_phenomatrix_list: Argument must be a PhenomatrixBase");
    }
#endif

    list<PhenomatrixBase> s;
    stack<Phenomatrix> p;

    // Allow different distance functions to be subbed in.
    double (*distance_function)(size_t, size_t, size_t, size_t);
};

#ifndef MATRIX_LIST_DEFINED
# define MATRIX_LIST_DEFINED
typedef std::list<PhenomatrixPair>          matrix_list;
#endif

#endif

