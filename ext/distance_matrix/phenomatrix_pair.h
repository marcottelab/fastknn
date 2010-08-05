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
#include <cmath>
#include <vector>
using std::vector;
using std::log;
using std::sqrt;
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
#include "similarity.h"
#include "phenomatrix.h"

#ifndef PHENOMATRIX_BASE_LIST_DEFINED
# define PHENOMATRIX_BASE_LIST_DEFINED
class Phenomatrix;
class PhenomatrixBase;
typedef list<PhenomatrixBase> phenomatrix_base_list;
#endif

class id_dist;
class id_dist_iter;
typedef std::priority_queue<id_dist_iter, std::vector<id_dist_iter>, std::greater<vector<id_dist_iter>::value_type> > proximity_queue;
typedef pair<sparse_document_vector, sparse_document_vector> sparse_document_vectors;

#ifndef MATRIX_LIST_DEFINED
# define MATRIX_LIST_DEFINED
class PhenomatrixPair;
typedef std::list<PhenomatrixPair>          matrix_list;
#endif

class PhenomatrixPair;

double hypergeometric(const PhenomatrixPair* const, uint, uint);
double manhattan(const PhenomatrixPair* const, uint, uint);
double euclidean(const PhenomatrixPair* const, uint, uint);
double cosine_similarity(const PhenomatrixPair* const, uint, uint);
double tanimoto_coefficient(const PhenomatrixPair* const, uint, uint);
double jaccard(const PhenomatrixPair* const, uint, uint);
double hellinger(const PhenomatrixPair* const, uint, uint);

// Handles creation of phenomatrices for pairs of species such that they have the
// correct number of rows.
//
// Also allows "temporary" removal of rows for the purpose of calculating
// distances and such (see push_mask and pop_mask).
class PhenomatrixPair {
public:
#ifdef RICE
    // Specialty constructor just for Ruby. Allows Ruby to pass in
    PhenomatrixPair(Rice::Object self, Rice::Object, Rice::Object, size_t);
#endif
    PhenomatrixPair(uint, PhenomatrixBase, size_t);

    PhenomatrixPair(uint id, uint given_id, size_t min_genes = 2)
    : s(create_phenomatrix_list(given_id, min_genes)),            // source species matrix
      p(create_phenomatrix_stack(id, given_id, min_genes)),       // predict species matrix
      distance_function(switch_distance_function("hypergeometric")),
      distance_threshold_(0.0)
    { }

    // Copy constructor
    PhenomatrixPair(const PhenomatrixPair& rhs)
    : s(rhs.s),
      p(rhs.p),
      distance_function(rhs.distance_function),
      distance_threshold_(rhs.distance_threshold_)
    { }

    ~PhenomatrixPair() { }
    
    uint id() const;
    uint predict_id() const;
    void push_mask(id_set mask_rows);
    bool pop_mask();


    // Return the distance, according to our distance function, between the two
    // columns. distance_threshold_ is probably only going to apply to distance
    // measures that utilize tf-idf (it's the idf threshold).
    double distance(uint j1, uint j2) const {
        return (*distance_function)( this, j1, j2 );
    }

    pair<size_t,size_t> observations_sizes(uint j1, uint j2) const;
    size_t total_column_count() const;
    size_t total_term_count(uint i) const;

    // Inverse Document Frequency for TF-IDF. Discretizes: threshold, by default,
    // is 0.0; but if you set it higher, everything that does not meet the threshold
    // will be treated as 0.
    double inverse_document_frequency(uint i) const {
        double idf = log( total_column_count() / (double)(total_term_count(i)) );
        if (idf < distance_threshold_) {
            //cerr << "phenomatrix_pair.h: inverse_document_frequency(i): distance_threshold_ = " << distance_threshold_ << "; returning 0" << endl;
            return 0.0;
        }
        else {
            //cerr << "phenomatrix_pair.h: inverse_document_frequency(i): distance_threshold_ = " << distance_threshold_ << "; returning idf = " << idf << endl;
            return idf;
        }
    }


    sparse_document_vector source_document_vector(uint) const;
    sparse_document_vector predict_document_vector(uint) const;
    sparse_document_vectors document_vectors(uint, uint) const;
    bool source_matrix_has_column(uint j) const;
    bool predict_matrix_has_column(uint j) const;

    // Find the single nearest neighbor (FIRST FOUND, not equivalence class)
    id_dist_iter nearest(uint j) const;


    // This is more of a helper for knearest. It needs to be set up properly in
    // DistanceMatrix.
    void knearest(proximity_queue& q, const uint& j, const size_t& k, double& kth_so_far, matrix_list::const_iterator this_iter) const;
    
    id_set intersection(uint j1, uint j2) const;
    size_t max_intersection_size() const;
    pair<size_t, size_t> min_observations_count() const;
    id_set observations(uint j) const;

    phenomatrix_base_list::const_reverse_iterator s_rbegin() const {
        return s.rbegin();
    }

    // Return the number of matrices on the stack.
    size_t size() const {
        return s.size();
    }

    // The lowest min_genes setting
    pair<size_t,size_t> min_genes() const;
    id_set row_ids() const;
    id_set predictable_column_ids() const;

#ifdef RICE
    Rice::Object get_distance_function() const {
        std::map<double(*)(const PhenomatrixPair* const, uint, uint), std::string> choices;
        choices[&hypergeometric]     = "hypergeometric";
        choices[&euclidean]          = "euclidean";
        choices[&manhattan]          = "manhattan";
        choices[&jaccard]            = "jaccard";
        choices[&hellinger]          = "hellinger";
        choices[&cosine_similarity]  = "cosine_similarity";
        choices[&tanimoto_coefficient] = "tanimoto_coefficient";

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

    void set_distance_threshold(float t) { distance_threshold_ = t; }
    float distance_threshold() const { return distance_threshold_; }
protected:

    // Calculate TF-IDF given some observations-size. Should only be used internally.
    double tf_idf_internal(uint i, uint j, size_t obs_size) const {
        return inverse_document_frequency(i) / (double)(obs_size);
    }

    // Returns a function pointer to a distance function based on a request made via
    // a string.
    static double (*switch_distance_function(const std::string& distance_measure))(const PhenomatrixPair* const, uint, uint) {
        std::map<std::string, double(*)(const PhenomatrixPair* const, uint, uint)> choices;
        choices["hypergeometric"]       = &hypergeometric;
        choices["euclidean"]            = &euclidean;
        choices["manhattan"]            = &manhattan;
        choices["jaccard"]              = &jaccard;
        choices["hellinger"]            = &hellinger;
        choices["cosine_similarity"]    = &cosine_similarity;
        choices["tanimoto_coefficient"] = &tanimoto_coefficient;

        return choices[distance_measure];
    }

    static stack<Phenomatrix> create_phenomatrix_stack(uint, uint, size_t);
    static list<PhenomatrixBase> create_phenomatrix_list(uint, size_t);
    static list<PhenomatrixBase> create_phenomatrix_list(const PhenomatrixBase&);
#ifdef RICE
    static stack<Phenomatrix> create_phenomatrix_stack(Rice::Object, uint, size_t);
    static list<PhenomatrixBase> create_phenomatrix_list(Rice::Object, size_t);
#endif

    list<PhenomatrixBase> s;
    stack<Phenomatrix> p;

    // Allow different distance functions to be subbed in.
    double (*distance_function)(const PhenomatrixPair* const, uint, uint);
    float distance_threshold_;
};


#ifndef MATRIX_LIST_DEFINED
# define MATRIX_LIST_DEFINED
typedef std::list<PhenomatrixPair>          matrix_list;
#endif

#endif

