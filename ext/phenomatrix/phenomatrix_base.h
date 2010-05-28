#ifndef PHENOMATRIX_BASE_H_
# define PHENOMATRIX_BASE_H_

typedef unsigned int uint;
#ifdef RICE
# include <rice/Object.hpp>
# include <rice/Data_Type.hpp>
# include <rice/Constructor.hpp>
#endif

#include <iostream>
//#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "connection.h"


#include <boost/unordered_map.hpp>

using std::cerr;
using std::cout;
using std::endl;
typedef boost::unordered_map<uint, std::set<uint> > omatrix;



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
    PhenomatrixBase(uint id, bool is_base_class = true)
    : id_(id), type_("Matrix"), obs(NULL)
    {
        base_construct(is_base_class);
    }


    // Copy constructor
    PhenomatrixBase(const PhenomatrixBase& rhs)
    : id_(rhs.id_),
      row_ids_(rhs.row_ids_),
      column_ids_(rhs.column_ids_),
      root_id_(rhs.root_id_),
      parent_id_(rhs.parent_id_),
      type_(rhs.type_),
      obs(new omatrix(*(rhs.obs)))
    {
#ifdef DEBUG_TRACE_COPY
        cerr << "phenomatrix_base.h: Copy constructor called! id = " << id_ << endl;
#endif
    }

    PhenomatrixBase(const PhenomatrixBase& rhs, const id_set& remove_rows)
    : id_(rhs.id_),
      row_ids_(rhs.row_ids_),
      column_ids_(rhs.column_ids_),
      root_id_(rhs.root_id_),
      parent_id_(rhs.parent_id_),
      type_(rhs.type_),
      obs(new omatrix(column_ids_.size()))
    {
#ifdef DEBUG_TRACE_COPY
        cerr << "phenomatrix_base.h: Copy constructor with mask called! id = " << id_ << endl;
#endif
        for (id_set::const_iterator c = column_ids_.begin(); c != column_ids_.end(); ++c) {
            // Copy the sets with certain things removed.
            set_difference((*(rhs.obs))[*c].begin(), (*(rhs.obs))[*c].end(),
                           remove_rows.begin(), remove_rows.end(),
                           std::insert_iterator<id_set>((*obs)[*c], (*obs)[*c].begin()));
        }
    }

    virtual ~PhenomatrixBase() {
        if (obs) delete obs;
    }

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
    size_t column_count() const { return obs->size(); }


    uint id() const { return id_; }

    // Get the ID of the earliest ancestor of this matrix
    uint root_id() const { return root_id_; }


    uint parent_id() const { return parent_id_; }

    // Get column
    id_set observations(uint j) const {
        omatrix::const_iterator jt = obs->find(j);
        if (jt == obs->end()) {
            cerr << "phenomatrix_base.h: observations(1): Requested non-existent column " << j << " on matrix " << id_ << endl;
            if (j == 0) throw;

#ifdef DEBUG_TRACE_DISTANCE
            // TEST CODE
            cerr << "Matrix is " << id() << endl;
            for (omatrix::const_iterator k = obs->begin(); k != obs->end(); ++k)
                cerr << k->first << ", " << std::flush;
            cerr << endl;
            // END TEST CODE
#endif
            string err = "phenomatrix_base.h: observations(1): Requested non-existent column " + lexical_cast<string>(j) + " on matrix " + lexical_cast<string>(id_);
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
        omatrix::const_iterator jt = obs->find(j);
        return (jt == obs->end()) ? false : (jt->second.find(i) != jt->second.end());
    }

    // Return the column identifiers
    id_set column_ids() const { return column_ids_; }
    // Return the row identifiers
    id_set row_ids() const { return row_ids_; }

    bool has_column(uint j) const {
#ifdef DEBUG_TRACE_DISTANCE
	cerr << "phenomatrix.h: has_column: on matrix " << id_ << ", requested col " << j << " and result will be " << (obs->find(j) != obs->end()) << endl;
#endif
        return (obs->find(j) != obs->end());
    }

    bool has_row(uint i) const {
        return (row_ids_.find(i) != row_ids_.end());
    }

    //bool has_observation(uint i, uint j) const {
    //    return ((*obs)[j].find(i));
    //}
    
protected:
    // Surrogate constructor, called by both default constructors, but not the
    // copy constructor.
    void base_construct(bool base = true) {
        // Get the information about the immediate parent and the root of the tree of matrices.
        std::pair<uint,uint> parent_and_root = fetch_parent_and_root_id();
        parent_id_ = parent_and_root.first;
        root_id_   = parent_and_root.second;

        if (base)
            inherit_construct();
    }

    // This will vary from base to child.
    void inherit_construct() {
        // Get matrix attributes
        row_ids_ = fetch_row_ids();
        column_ids_ = fetch_column_ids();

        // Create the observation matrix
        obs = new omatrix(column_ids_.size());

        load_matrix();
    }
    

    void add_observation(uint i, uint j) {
        (*obs)[j].insert(i);
    }

    void remove_observation(uint i, uint j) {
        (*obs)[j].erase(i);
    }


    // There are other ways to arrange these two functions, but only this one appears
    // to work with Rice.
    string fetch_matrix_type(uint matrix_id) {
        return Connection::instance().fetch_type("matrices", matrix_id);
    }
    string fetch_matrix_type() {
        return Connection::instance().fetch_type("matrices", id_);
    }


    virtual string load_matrix_sql(uint matrix_id) const {
#ifdef DEBUG_TRACE_INHERITED_CONSTRUCTION
        cerr << "phenomatrix_base.h: load_matrix_sql" << endl;
#endif
        return "SELECT DISTINCT i,j FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + " AND type = 'Cell';";
    }
    virtual string row_count_sql(uint matrix_id) const {
        return "SELECT COUNT(DISTINCT i) FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + ";";
    }
    virtual string row_ids_sql(uint matrix_id) const {
        return "SELECT DISTINCT i FROM entries WHERE matrix_id = " + lexical_cast<string>(matrix_id) + " ORDER BY i;";
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

        for (result_t::const_iterator rt = r.begin(); rt != r.end(); ++rt) {
            uint i = 0, j = 0;
            (*rt)[0].to(i);
            (*rt)[1].to(j);

            add_observation(i, j);
        }

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
        return Connection::instance().fetch_count("SELECT COUNT(DISTINCT j) FROM entries WHERE matrix_id = " + lexical_cast<string>(of_id) + ";");
    }

    id_set fetch_column_ids(uint of_id) const {
        return Connection::instance().fetch_id_set("SELECT DISTINCT j FROM entries WHERE matrix_id = " + lexical_cast<string>(of_id) + " ORDER BY j;");
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


    uint id_;

    id_set row_ids_;
    id_set column_ids_;
    uint root_id_;
    uint parent_id_;
    string type_;

    omatrix* obs; // set of observations indexed by column
};

#endif // PHENOMATRIX_H_
