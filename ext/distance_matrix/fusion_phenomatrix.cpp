#include "fusion_phenomatrix.h"

id_set extract_matrix_ids(const matrix_list& matrix_pairs) {
    id_set ret;
    BOOST_FOREACH( PhenomatrixPair matrix_pair, matrix_pairs ) {
        ret.insert(matrix_pair.id());
    }
    return ret;
}


id_set extract_row_ids(const matrix_list& matrix_pairs) {

}

id_set extract_column_ids(const matrix_list& matrix_pairs) {

}


FusionPhenomatrix::FusionPhenomatrix(uint id, const matrix_list& given_matrices)
: PhenomatrixBase(id, false),
  given_ids_(extract_matrix_ids(given_matrices)) {
    row_ids_    = extract_row_ids(given_matrices);
    column_ids_ = extract_column_ids(given_matrices);
}