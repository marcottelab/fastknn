#ifndef PHENOMATRIX_H_
# define PHENOMATRIX_H_

typedef unsigned int uint;
#undef ALLOC // Needed for annoying bug (not in my code)
#include <boost/numeric/ublas/vector_sparse.hpp>
typedef boost::numeric::ublas::compressed_vector<double> sparse_document_vector;
#undef ALLOC // Needed for annoying bug (not in my code)

#ifdef RICE
# include <rice/Object.hpp>
# include <rice/Data_Type.hpp>
# include <rice/Constructor.hpp>
#endif
#include <map>
#include <queue>
#include <iostream>
#include <cmath>
#include <string>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>
using std::log;
using std::map;
using std::queue;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using boost::lexical_cast;

#include "connection.h"
#include "phenomatrix_pair.h"

typedef boost::unordered_map<uint, std::set<uint> > omatrix;
typedef boost::unordered_map<uint, size_t> gene_counter;

class PhenomatrixPair;

// Stores a single species phenomatrix without concern for the rows in the source
// species.
//
// Remember, the source species isn't going to have as many rows as the predict
// species. Phenomatrix, which inherits from PhenomatrixBase, will take into
// account both species when producing total row counts, etc.
class PhenomatrixBase {
public:
    // This constructor is typically called by DistanceMatrix, so different matrices
    // can share the same database connection.
    //
    // At some point, it should be revised to accept a Ruby on Rails connection.
    PhenomatrixBase(uint id, bool is_base_class = true, size_t min_genes = 2)
    : id_(id), min_genes_(min_genes), child_ids_(fetch_child_ids()), type_("Matrix")
    {
        base_construct(is_base_class, min_genes);
    }

    PhenomatrixBase(const PhenomatrixBase& rhs);
    PhenomatrixBase(const PhenomatrixBase& rhs, const id_set& remove_rows);
    virtual ~PhenomatrixBase() { }

#ifdef RICE
    Rice::Object rb_root_id() const {
        uint rid = root_id();
        return (rid == 0) ? Rice::Object() : Rice::Object(to_ruby<uint>(rid));
    }
    Rice::Object rb_parent_id() const {
        uint pid = parent_id();
        return (pid == 0) ? Rice::Object() : Rice::Object(to_ruby<uint>(pid));
    }
#endif

    // Get the number of unique rows in the matrix. Naive to node/tree/leaf status.
    size_t row_count() const { return row_ids_.size(); }

    // Returns the number of columns in this matrix tree, which may not be the
    // same as the number of columns in the obs matrix (but should).
    //size_t tree_column_count() const { return column_count_; }

    // Returns the number of columns in the fully-loaded matrix -- may be different
    // than expected if this is a mask.
    // This may also used as the number of documents for TF-IDF.
    size_t column_count() const { return obs.size(); }

    // Returns the number of times some gene is found in the matrix. Does not
    // take into account any masks that might have been pushed, since this is
    // just for TF-IDF.
    size_t term_count(uint i) const {
        gene_counter::const_iterator it = gcount.find(i);
        if (it == gcount.end())  return 0;
        else                     return it->second;
    }

    // Inverse Document Frequency for TF-IDF. Discretizes: threshold, by default,
    // is 0.0; but if you set it higher, everything that does not meet the threshold
    // will be treated as 0.
/*    double inverse_document_frequency(uint i, float idf_threshold) const {
        double idf = log( column_count() / (double)(term_count(i)) );
        if (idf < idf_threshold) return 0.0;
        else                     return idf;
    } */

    // Return a TF-IDF sparse vector for some column (phenotype).
    sparse_document_vector document_vector(uint, const PhenomatrixPair* const) const;

    uint id() const { return id_; }

    // Get the ID of the earliest ancestor of this matrix
    uint root_id() const { return root_id_; }


    uint parent_id() const { return parent_id_; }

    // Return the min genes setting -- not the same as min_observations_count,
    // which returns the actual minimum number of genes in any column.
    size_t min_genes() const { return min_genes_; }

    size_t min_observations_count() const {
        size_t min = UINT_MAX;
        for (omatrix::const_iterator it = obs.begin(); it != obs.end(); ++it)
            if (it->second.size() < min) min = it->second.size();
        return min;
    }

    // Get column
    id_set observations(uint j) const {
        omatrix::const_iterator jt = obs.find(j);
        if (jt == obs.end()) {

#ifdef DEBUG_TRACE_DISTANCE
            // TEST CODE
            cerr << "Matrix is " << id() << endl;
            for (omatrix::const_iterator k = obs.begin(); k != obs.end(); ++k)
                cerr << k->first << ", " << std::flush;
            cerr << endl;
            // END TEST CODE
#endif
            string err = "phenomatrix_base.h: observations: Requested non-existent column " + lexical_cast<string>(j) + " on matrix " + lexical_cast<string>(id_);
#ifdef RICE
            throw Rice::Exception(rb_eArgError, err.c_str());
#else
            cerr << err << endl;
            throw;
#endif
        }
        return jt->second;
    }

    // Get column size
    size_t observations_size(uint j) const { return observations(j).size(); }

    string type() const { return type_; }

    // Determine whether there's an observation at gene i, phene j
    bool operator()(uint i, uint j) const {
        omatrix::const_iterator jt = obs.find(j);
        return (jt == obs.end()) ? false : (jt->second.find(i) != jt->second.end());
    }

    // Return the column identifiers
    id_set column_ids() const { return column_ids_; }
    // Return the row identifiers
    id_set row_ids() const { return row_ids_; }

    bool has_column(uint j) const {
#ifdef DEBUG_TRACE_DISTANCE
	cerr << "phenomatrix.h: has_column: on matrix " << id_ << ", requested col " << j << " and result will be " << (obs.find(j) != obs.end()) << endl;
#endif
        return (column_ids_.find(j) != column_ids_.end());
    }

    bool has_row(uint i) const {
        return (row_ids_.find(i) != row_ids_.end());
    }

    //bool has_observation(uint i, uint j) const {
    //    return ((*obs)[j].find(i));
    //}


    // This function returns the ids for the matrix's children. To get those
    // IDs requires a database fetch -- it is not pre-cached on load.
    id_set child_ids() const {
        id_set cids;
        for (map<uint,string>::const_iterator i = child_ids_.begin(); i != child_ids_.end(); ++i)
            cids.insert(i->first);
        return cids;
    }

    std::map<uint, id_set> child_row_ids() const {
        return fetch_child_row_ids();
    }

protected:

    // Surrogate constructor, called by both default constructors, but not the
    // copy constructor.
    void base_construct(bool base, size_t min_genes) {
        // Get the information about the immediate parent and the root of the tree of matrices.
        std::pair<uint,uint> parent_and_root = fetch_parent_and_root_id();
        parent_id_ = parent_and_root.first;
        root_id_   = parent_and_root.second;

        if (base)
            inherit_construct(min_genes);
    }

    // This will vary from base to child.
    void inherit_construct(size_t min_genes) {
        // Get matrix attributes
        row_ids_ = fetch_row_ids();
        // column_ids_ = fetch_column_ids(); // Now set by enforce_min_genes

        // Create the observation matrix
        obs = omatrix(fetch_column_count());

        // Create the index giving the number of phenotypes where a gene appears
        gcount = gene_counter(row_ids_.size());

        load_matrix();

        // Set column_ids_ and make sure columns have at least min_genes observations.
        enforce_min_genes(min_genes);
    }


    // This will not vary at all -- removes columns from the matrix that have
    // below min_genes items in them.
    void enforce_min_genes(size_t min_genes) {
        // ensure that each column has at least min_genes in it, and create
        // column_ids_ based on that.
        queue<omatrix::const_iterator> erase_q;

        for (omatrix::iterator it = obs.begin(); it != obs.end(); ++it) {
            if (it->second.size() >= min_genes) {

                //cerr << "j=" << it->first << "\t min_genes = " << min_genes << "\tit->second.size() = " << it->second.size() << endl;
                column_ids_.insert(it->first);

                // Walk through and increment gcount for each gene found.
                for (std::set<uint>::const_iterator gt = it->second.begin(); gt != it->second.end(); ++gt) {
                    if (gcount.find(*gt) == gcount.end())   gcount[*gt] = 1;
                    else                                    gcount[*gt]++;
                }

            } else erase_q.push(it);
        }

        // Remove columns from obs that are below the threshold
        while (erase_q.size() > 0) {
            obs.erase(erase_q.front());
            erase_q.pop();
        }

        min_genes_ = min_genes;
    }


    void add_observation(uint i, uint j) {
        obs[j].insert(i);
    }

    void add_observations(id_set i_set, uint j) {
        obs[j] = i_set;
    }

    void remove_observation(uint i, uint j) {
        obs[j].erase(i);
    }


    // There are other ways to arrange these two functions, but only this one appears
    // to work with Rice.
    string fetch_matrix_type(uint matrix_id) const {
        return Connection::instance().fetch_type("matrices", matrix_id);
    }

    string fetch_matrix_type() const {
        try { // There may be a better place for this exception. We shall see.
            return Connection::instance().fetch_type("matrices", id_);
        } catch(...) {
            string err = "It appears that the given matrix id (" + lexical_cast<string>(id_) + ") is not in the database.";
#ifdef RICE
            throw Rice::Exception(rb_eArgError, err.c_str());
#else
            cerr << err << endl;
            throw;
#endif
        }
    }

    map<uint,string> fetch_child_ids(uint matrix_id) const {
        return Connection::instance().fetch_map<uint,string>(child_ids_sql(matrix_id));
    }

    map<uint,string> fetch_child_ids() const {
        return Connection::instance().fetch_map<uint,string>(child_ids_sql(id_));
    }

    // Same as fetch_child_row_ids, but treats the child as a node instead of a leaf.
    id_set fetch_node_as_leaf_row_ids(uint matrix_id, uint child_matrix_id) const {
        return Connection::instance().fetch_id_set(node_as_leaf_row_ids_sql(matrix_id, child_matrix_id));
    }
    id_set fetch_leaf_row_ids(uint of_id) const {
        return Connection::instance().fetch_id_set(leaf_row_ids_sql(of_id));
    }

    std::map<uint, id_set> fetch_child_row_ids() const {
        std::map<uint,id_set> child_to_row_ids;
        for (map<uint,string>::const_iterator ct = child_ids_.begin(); ct != child_ids_.end(); ++ct) {
            if (ct->second == "LeafMatrix")
                child_to_row_ids[ct->first] = fetch_leaf_row_ids(ct->first);
            else
                child_to_row_ids[ct->first] = fetch_node_as_leaf_row_ids(id_, ct->first);
        }
        return child_to_row_ids;
    }


    virtual string load_matrix_sql(uint matrix_id) const {
#ifdef DEBUG_TRACE_INHERITED_CONSTRUCTION
        cerr << "phenomatrix_base.h: load_matrix_sql" << endl;
#endif
        return "SELECT DISTINCT i,j FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + " AND type = 'Cell' ORDER BY j, i;";
    }
    virtual string row_count_sql(uint matrix_id) const {
        return "SELECT COUNT(DISTINCT i) FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + ";";
    }
    virtual string row_ids_sql(uint matrix_id) const {
        return leaf_row_ids_sql(matrix_id);
    }
    virtual string column_count_sql(uint matrix_id) const {
        return "SELECT COUNT(DISTINCT j) FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + ";";
    }
    virtual string column_ids_sql(uint matrix_id) const {
        return "SELECT DISTINCT j FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + " ORDER BY j;";
    }

    virtual string leaf_row_ids_sql(uint matrix_id) const {
        return "SELECT DISTINCT i FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + " ORDER BY i;";
    }

    // THESE ARE NOT DESIGNED TO ACCOUNT FOR source matrices in the same way as row_ids_sql
    // (which is overridden in phenomatrix.h)
    virtual string node_as_leaf_row_ids_sql(uint matrix_id, uint node_matrix_id) const {
        return "SELECT DISTINCT e1.i FROM entries e1 WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) +
               " EXCEPT SELECT e2.i FROM entries e2 WHERE e2.matrix_id = " + lexical_cast<string>(node_matrix_id) +
               " ORDER BY i;";
    }
    virtual string child_ids_sql(uint matrix_id) const {
        return "SELECT DISTINCT id, type FROM matrices WHERE parent_id = " + lexical_cast<string>(matrix_id) + " ORDER BY id;";
    }


    // Force a matrix to be loaded as a Node
    void load_matrix(uint matrix_id) {
        // Find out whether this is a node or leaf matrix.
        // type_ = fetch_type(matrix_id);

        //work_t w(*c);
        work_t* w = Connection::instance().work();
        result_t r = w->exec( load_matrix_sql(matrix_id) );

        if (r.size() == 0) {
            string err = "phenomatrix_base.h: load_matrix(1): Could not load matrix.";
#ifdef RICE
            throw Rice::Exception(rb_eArgError, err.c_str());
#else
            cerr << err << endl;
            throw;
#endif
        }


        id_set set_of_i_for_this_j;
        result_t::const_iterator rt = r.begin();
        if (rt == r.end()) {
            string err = "phenomatrix_base.h: load_matrix: Empty matrix!";
#ifdef RICE
            throw Rice::Exception(rb_eArgError, err.c_str());
#else
            cerr << err << endl;
            throw;
#endif
        }

        uint this_j = 0;
        (*rt)[0].to(this_j);
        id_set::iterator insert_hint;

        while (rt != r.end()) {
            uint i = 0, j = 0;
            (*rt)[0].to(i);
            (*rt)[1].to(j);

            // Any time we get to a
            if (this_j != j) {
                add_observations(set_of_i_for_this_j, this_j);
                set_of_i_for_this_j.clear();

                // Prepare to insert the next j's i-values
                insert_hint = set_of_i_for_this_j.begin();
                this_j = j;
            }
            insert_hint = set_of_i_for_this_j.insert(insert_hint, i);
            ++rt;
        }

        if (set_of_i_for_this_j.size() > 0)
            add_observations(set_of_i_for_this_j, this_j);

        delete w;
    }



    // Force a mask on top of a Node
    void mask_load_matrix() {
        work_t* w = Connection::instance().work();
        result_t r = w->exec( this->load_matrix_sql(id_) );

        for (result_t::const_iterator rt = r.begin(); rt != r.end(); ++rt) {
            uint i = 0, j = 0;
            (*rt)[0].to(i);
            (*rt)[1].to(j);

            remove_observation(i, j);
        }
        delete w;
    }


    // Read all entries from the matrix and store them in obs
    void load_matrix() {
        // If this matrix is a leaf, we need to load the parent first
        type_ = fetch_matrix_type();

        if (type_ == "LeafMatrix") { // This is not the most efficient way of doing things, but good enough.
            load_matrix(parent_id()); // Add cells
            mask_load_matrix();     // Then remove the masked ones
        } else {
            if (type_ == "") type_ = "Matrix";
            load_matrix(id_);
        }
    }


    id_set fetch_row_ids(uint of_id) const {
        return Connection::instance().fetch_id_set(row_ids_sql(of_id));
    }
    id_set fetch_row_ids() const { return fetch_row_ids(id_); }

    // Determine the number of unique rows in the matrix. Naive to node/tree/leaf status.
    size_t fetch_row_count(uint of_id) const {
        return Connection::instance().fetch_count(row_count_sql(of_id));
    }
    size_t fetch_row_count() const { return fetch_row_count(id_); } // default arg

    size_t fetch_column_count(uint of_id) const {
        return Connection::instance().fetch_count(column_count_sql(of_id));
    }

    id_set fetch_column_ids(uint of_id) const {
        return Connection::instance().fetch_id_set(column_ids_sql(of_id));
    }
    id_set fetch_column_ids() const { return fetch_column_ids(id_); }

    size_t fetch_column_count() const { // default arg
        uint id = id_;
        if (parent_id_ > 0) {
            if (root_id_ > 0) id = root_id_;
            else              id = parent_id_;
        }

        return fetch_column_count(id);
    }

    // Determine the number of rows in the matrix tree.
    size_t fetch_root_row_count() const {
        uint rid = root_id();
        if (rid == 0 || rid == id_)
            return row_ids_.size();
        else
            return Connection::instance().fetch_count("SELECT COUNT(DISTINCT i) FROM entries WHERE matrix_id = " + lexical_cast<string>(rid) + ";");
    }

    std::pair<uint,uint> fetch_parent_and_root_id() const {
        uint par;
        uint x = par = parent_id(id_);
        uint prev_x = id_;
        while (x != 0) {
            prev_x = x;
            x = parent_id(x);
        }
        return std::make_pair<uint,uint>(par,prev_x);
    }


    // Get the ID of the parent of some other matrix
    uint parent_id(uint of_id) const {
        return Connection::instance().fetch_id("SELECT DISTINCT parent_id FROM matrices WHERE id = " + lexical_cast<string>(of_id) + ";");
    }


    // COPY CONSTRUCTION HELPERS
    static omatrix copy_construct_omatrix(const PhenomatrixBase& rhs);
    static omatrix copy_construct_omatrix(const PhenomatrixBase& rhs, const id_set& remove_rows);

    static id_set remove_from_column(const id_set& rhs_col, const id_set& remove_rows) {
        id_set lhs_col;
        // Copy the sets with certain things removed.
        set_difference(rhs_col.begin(), rhs_col.end(),
                       remove_rows.begin(), remove_rows.end(),
                       std::insert_iterator<id_set>(lhs_col, lhs_col.begin()));
        return lhs_col;
    }

    static id_set copy_construct_row_ids(const PhenomatrixBase& rhs, const id_set& remove_rows);


    uint id_;
    uint min_genes_;

    id_set row_ids_;
    id_set column_ids_;
    uint root_id_;
    uint parent_id_;
    map<uint,string> child_ids_;
    string type_;

    omatrix obs;            // set of observations indexed by column

    gene_counter gcount;    // Count of observations indexed by row (for tf-idf)
                            // Note that gcount does not change when the matrix is masked!
};



// This derivative class of PhenomatrixBase will typically represent the prediction
// species in a pair of species. It can be used for either! Just set the two ids
// to the same value.
//
// As described in PhenomatrixBase's comment section, a prediction matrix will
// have more rows than a source matrix. This Phenomatrix takes that into account
// when loading the prediction matrix and only loads the part which is relevant
// to the source matrix.
class Phenomatrix : public PhenomatrixBase {
public:
    Phenomatrix(uint, uint, size_t);
    Phenomatrix(const Phenomatrix&);
    Phenomatrix(const Phenomatrix&, const id_set&);
    ~Phenomatrix() { }

    // Get the source matrix id
    uint source_id() const { return given_id_; }
protected:
    // OVERRIDE SQL FROM PHENOMATRIX_BASE.
    string load_matrix_sql(uint matrix_id) const {
        return string("SELECT DISTINCT e1.i,e1.j FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + " \n" +
               " AND  e1.type = 'Cell' ORDER BY e1.j, e1.i;";
    }

    string row_count_sql(uint matrix_id) const {
        return string("SELECT COUNT(DISTINCT e1.i) FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + ";";
    }

    string row_ids_sql(uint matrix_id) const {
        return string("SELECT DISTINCT e1.i FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + " ORDER BY e1.i;";
    }

    string column_ids_sql(uint matrix_id) const {
        return string("SELECT DISTINCT e1.j FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + " ORDER BY e1.j;";
    }

    string column_count_sql(uint matrix_id) const {
        return string("SELECT COUNT(DISTINCT e1.j) FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + ";";
    }

    uint given_id_;
};

#endif
