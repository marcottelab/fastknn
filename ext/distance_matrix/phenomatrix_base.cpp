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

    sparse_document_vector v(row_ids_.size(), j_obs.size());
    for (id_set::const_iterator it = j_obs.begin(); it != j_obs.end(); ++it) {
        v[*it] = idf_owner->inverse_document_frequency(*it) / (double)(j_obs.size());
    }

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