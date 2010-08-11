//#include "phenomatrix_base.h"

// Copy constructor
PhenomatrixBase::PhenomatrixBase(const PhenomatrixBase& rhs)
: id_(rhs.id_),
  min_genes_(rhs.min_genes_),
  row_ids_(rhs.row_ids_),
  column_ids_(rhs.column_ids_),
  root_id_(rhs.root_id_),
  parent_id_(rhs.parent_id_),
  child_ids_(rhs.child_ids_),
  type_(rhs.type_),
  obs(copy_construct_omatrix(rhs)),
  gcount(rhs.gcount)
{ }

PhenomatrixBase::PhenomatrixBase(const PhenomatrixBase& rhs, const id_set& remove_rows)
: id_(rhs.id_),
  min_genes_(rhs.min_genes_),
  row_ids_(copy_construct_row_ids(rhs, remove_rows)),
  column_ids_(rhs.column_ids_),
  root_id_(rhs.root_id_),
  parent_id_(rhs.parent_id_),
  child_ids_(rhs.child_ids_), // not safe to use!
  type_(rhs.type_),
  obs(copy_construct_omatrix(rhs, remove_rows)),
  gcount(rhs.gcount)
{ }


// Return a TF-IDF sparse vector for some column (phenotype).
sparse_document_vector PhenomatrixBase::document_vector(uint j, const PhenomatrixPair* const idf_owner) const {
    id_set j_obs = observations(j);

    // This is not the cleanest container type to use here. Initially, I had set
    // the size to row_ids_.size(), but that was problematic because certain genes (esp. for At)
    // had IDs way above that size. These were getting lost! Now, it sets size
    // as the maximum row id.
    // compressed_vector is a component of compressed_matrix, which uses "Compressed Row Storage"
    // from NetLib: http://www.netlib.org/linalg/html_templates/node91.html
    // Based on that piece of information, it looks as if empty rows do not matter
    // at all. As such, it should be perfectly safe to set an insanely high size,
    // though I wonder why this type needs one at all.
    // Oops. What is the type on unbounded_array? I have a question pending on stackoverflow,
    // but for now I'll leave this as is.
    sparse_document_vector v(*(row_ids_.rbegin()) + 1, j_obs.size());
    for (id_set::const_iterator it = j_obs.begin(); it != j_obs.end(); ++it) {
        cerr << "phenomatrix_base.cpp: document_vector: IDF for " << *it << "\t" << idf_owner->inverse_document_frequency(*it) << "\tdenom: " << j_obs.size() << endl;
        double d = idf_owner->inverse_document_frequency(*it) / (double)(j_obs.size());
        v.insert_element(*it, d); // faster than v[*it], I speculate based on: http://www.guwi17.de/ublas/matrix_sparse_usage.html
        // cerr << "v: " << *it << "\t" << const_castv[*it] << endl; // REMOVE ME
    }
/*    cerr << "Next:" << endl;

    for (sparse_document_vector::const_iterator it = v.begin(); it != v.end(); ++it) {
        cerr << "v: " << it.index() << "\t" << *it << endl; // REMOVE ME
    }
    cerr << "j_obs size = " << j_obs.size() << endl; */

    return v;
}


// COPY CONSTRUCTION HELPERS
omatrix PhenomatrixBase::copy_construct_omatrix(const PhenomatrixBase& rhs) {
    return omatrix(rhs.obs);
}

omatrix PhenomatrixBase::copy_construct_omatrix(const PhenomatrixBase& rhs, const id_set& remove_rows) {
    omatrix lhs(rhs.column_ids_.size());

    // Remove rows from each column
    for (id_set::const_iterator c = rhs.column_ids_.begin(); c != rhs.column_ids_.end(); ++c) {
        id_set rhs_col = rhs.obs.find(*c)->second;
        lhs[*c] = remove_from_column(rhs_col, remove_rows);
    }
    return lhs;
}

id_set PhenomatrixBase::copy_construct_row_ids(const PhenomatrixBase& rhs, const id_set& remove_rows) {
    id_set lhs_row_ids;
    // Remove rows from row_ids list
    set_difference(rhs.row_ids_.begin(), rhs.row_ids_.end(),
                   remove_rows.begin(), remove_rows.end(),
                   std::insert_iterator<id_set>(lhs_row_ids, lhs_row_ids.begin()));
    return lhs_row_ids;
}