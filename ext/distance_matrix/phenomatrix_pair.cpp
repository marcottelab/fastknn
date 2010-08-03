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


// Cosine similarity between the two matrices in pair p, using TF-IDF, and
// treating both as a single corpus.
double cosine_similarity(const PhenomatrixPair* const p, uint j1, uint j2, float threshold) {
    double dot_product_accum = 0.0;
    double p_magnitude_accum = 0.0;
    double s_magnitude_accum = 0.0;

    id_set row_ids = p->row_ids();
    for (id_set::const_iterator it = row_ids.begin(); it != row_ids.end(); ++it) {

        double p_component = p->tf_idf(*it, j1), s_component = p->tf_idf(*it, j2);

        dot_product_accum += p_component * s_component;
        p_magnitude_accum += p_component * p_component;
        s_magnitude_accum += s_component * s_component;
    }

    double sim = dot_product_accum / sqrt(p_magnitude_accum * s_magnitude_accum);
    if (sim > 1.0) return 0.0;
    else return 1.0 - sim;
}

// Tanimoto coefficient between the two matrices in pair p, using TF-IDF, and
// treating both as a single corpus.
double tanimoto_coefficient(const PhenomatrixPair* const p, uint j1, uint j2, float threshold) {
    double dot_product_accum = 0.0;
    double p_magnitude_accum = 0.0;
    double s_magnitude_accum = 0.0;

    for (id_set::const_iterator it = p->row_ids().begin(); it != p->row_ids().end(); ++it) {

        double p_component = p->tf_idf(*it, j1), s_component = p->tf_idf(*it, j2);

        dot_product_accum += p_component * s_component;
        p_magnitude_accum += p_component * p_component;
        s_magnitude_accum += s_component * s_component;
    }

    double sim = dot_product_accum / (p_magnitude_accum + s_magnitude_accum - dot_product_accum);
    if (sim > 1.0) return 0.0;
    else return 1.0 - sim;
}


double jaccard(const PhenomatrixPair* const p, uint j1, uint j2, float threshold) {
    pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;

    return jaccard(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );
}


double hellinger(const PhenomatrixPair* const p, uint j1, uint j2, float threshold) {
    pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;

    return hellinger(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );
}

// Exposes manhattan(m,n,k,N) and PhenomatrixPair to eachother.
double manhattan(const PhenomatrixPair* const p, uint j1, uint j2, float threshold) {
    pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;

    return manhattan(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );
}

// Exposes euclidean(m,n,k,N) and PhenomatrixPair to eachother.
double euclidean(const PhenomatrixPair* const p, uint j1, uint j2, float threshold) {
    pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;

    return euclidean(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );
}

// Exposes hypergeometric(m,n,k,N) and PhenomatrixPair to eachother.
// Each distance function should have something like this -- and it's important
// that it have exactly these arguments.
double hypergeometric(const PhenomatrixPair* const p, uint j1, uint j2, float threshold) {
    std::pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;

    return hypergeometric(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );
}