#include "fusion_phenomatrix.h"

id_set extract_matrix_ids(const matrix_list& matrix_pairs) {
    id_set ret;
    BOOST_FOREACH( const PhenomatrixPair& matrix_pair, matrix_pairs ) {
        ret.insert(matrix_pair.id());
    }
    return ret;
}

id_set intersect_ids(id_set a, id_set b) {
    id_set ret;
    set_intersection(a.begin(), a.end(),
                     b.begin(), b.end(),
                     insert_iterator<id_set>(ret, ret.begin()));
    return ret;
}

id_set union_ids(id_set a, id_set b) {
    id_set ret;
    set_union(a.begin(), a.end(),
              b.begin(), b.end(),
              insert_iterator<id_set>(ret, ret.begin()));
    return ret;
}

id_set extract_row_ids(const matrix_list& matrix_pairs) {
    id_set ret;
    BOOST_FOREACH( const PhenomatrixPair& matrix_pair, matrix_pairs ) {
        ret = union_ids(matrix_pair.row_ids(), ret);
    }
    return ret;
}

id_set extract_column_ids(const matrix_list& matrix_pairs) {
    id_set ret;
    BOOST_FOREACH( const PhenomatrixPair& matrix_pair, matrix_pairs ) {
        ret = union_ids(matrix_pair.predictable_column_ids(), ret);
    }
    return ret;
}


FusionPhenomatrix::FusionPhenomatrix(uint id, const matrix_list& given_matrices)
: PhenomatrixBase(id, false),
  given_ids_(extract_matrix_ids(given_matrices)) {
    row_ids_    = extract_row_ids(given_matrices);
    column_ids_ = extract_column_ids(given_matrices);
}