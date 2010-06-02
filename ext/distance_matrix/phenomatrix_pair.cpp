#include "phenomatrix_base.h"
#include "phenomatrix.h"
#include "phenomatrix_pair.h"


// Find the single nearest neighbor (FIRST FOUND, not equivalence class)
id_dist_iter PhenomatrixPair::nearest(uint j) const {

    double min_dist = 100;
    uint min_dist_id = 0;

    id_set s_columns = s.back().column_ids();

    for (id_set::const_iterator k = s_columns.begin(); k != s_columns.end(); ++k) {
        if (j == *k) continue; // Don't count it when the columns are the same
        double d_jk = distance(j, *k);
        if (d_jk < min_dist) {
            min_dist = d_jk;
            min_dist_id = *k;
        }
    }

    return id_dist_iter(min_dist_id, min_dist, id());
}


// This is more of a helper for knearest. It needs to be set up properly in
// DistanceMatrix.
void PhenomatrixPair::knearest(proximity_queue& q, const uint& j, const size_t& k, double& kth_so_far, matrix_list::const_iterator this_iter) const {
    const id_set& columns = s.back().column_ids();

    uint current_k = 0;

    // Only sort all if this is the first matrix added (e.g., if kth_so_far
    // still very high). After that, we can limit to those within kth_so_far.
    for (id_set::const_iterator st = columns.begin(); st != columns.end(); ++st) {
        if (j == *st) continue;
        double d_jk = distance(j, *st);
        if (d_jk <= kth_so_far) { // don't add distances that are already outside the boundary (kth_so_far)
            q.push( id_dist_iter(*st, d_jk, this_iter) );

            // We can now change the boundary to be whatever is the kth value,
            // which will save us time on the next matrix.
            ++current_k;
            if (current_k == k) kth_so_far = d_jk;
        }
    }
}

