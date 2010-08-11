#include "phenomatrix.h"
#include "phenomatrix_pair.h"


#ifdef RICE
// Specialty constructor just for Ruby. Allows Ruby to pass in
PhenomatrixPair::PhenomatrixPair(Rice::Object self, Rice::Object predict_or_id, Rice::Object source_or_id, size_t min_genes)
: s(create_phenomatrix_list(source_or_id, min_genes)),
  p(create_phenomatrix_stack(predict_or_id, s.back().id(), min_genes)),
  distance_function(switch_distance_function("hypergeometric")),
  min_idf_(0.0)
{ }
#endif

PhenomatrixPair::PhenomatrixPair(uint id, PhenomatrixBase given_phenomatrix, size_t min_genes = 2)
: s(create_phenomatrix_list(given_phenomatrix)),
  p(create_phenomatrix_stack(id, given_phenomatrix.id(), min_genes)),
  distance_function(switch_distance_function("hypergeometric")),
  min_idf_(0.0)
{ }


// Return the ID of the source matrix
uint PhenomatrixPair::id() const {
    return s.back().id();
}

uint PhenomatrixPair::predict_id() const {
    return p.top().id();
}


// Removes rows from the matrices on which we're calculating distances.
void PhenomatrixPair::push_mask(id_set mask_rows) {
    p.push(      Phenomatrix(     p.top() , mask_rows ) );
    s.push_back( PhenomatrixBase( s.back(), mask_rows ) );
}


// Restores removed rows from the matrices on which we're calculating distances.
bool PhenomatrixPair::pop_mask() {
    if (p.size() == 1) return false;
    p.pop();
    s.pop_back();
    return true;
}


pair<size_t,size_t> PhenomatrixPair::observations_sizes(uint j1, uint j2) const {
    return make_pair<size_t,size_t>(p.top().observations_size(j1), s.back().observations_size(j2));
}

// Document count for TF-IDF when we consider the species pair to be the
// corpus (rather than each species as its own corpus).
size_t PhenomatrixPair::total_column_count() const {
    return s.back().column_count() + p.top().column_count();
}

// Total term count for a given term, for TF-IDF, when we consider the species
// pair to be the corpus (rather than each species as its own corpus).
size_t PhenomatrixPair::total_term_count(uint i) const {
    return s.back().term_count(i) + p.top().term_count(i);
}


sparse_document_vector PhenomatrixPair::source_document_vector(uint j) const {
    return s.back().document_vector(j, this);
}
sparse_document_vector PhenomatrixPair::predict_document_vector(uint j) const {
    return p.top().document_vector(j, this);
}

sparse_document_vectors PhenomatrixPair::document_vectors(uint j1, uint j2) const {
    sparse_document_vector v1(p.top().document_vector(j1, this));
    sparse_document_vector v2(s.back().document_vector(j2, this));

    return make_pair(v1, v2);
}


bool PhenomatrixPair::source_matrix_has_column(uint j) const {
    return s.back().has_column(j);
}

bool PhenomatrixPair::predict_matrix_has_column(uint j) const {
    return p.top().has_column(j);
}


id_set PhenomatrixPair::intersection(uint j1, uint j2) const {
    const id_set& s1 = p.top().observations(j1);
    const id_set& s2 = s.back().observations(j2);

    id_set ret;
    set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                      std::insert_iterator<id_set>(ret,ret.begin()));
    return ret;
}


size_t PhenomatrixPair::max_intersection_size() const {
    return s.back().row_count();
}

pair<size_t, size_t> PhenomatrixPair::min_observations_count() const {
    return make_pair<size_t,size_t>(p.top().min_observations_count(), s.back().min_observations_count());
}

id_set PhenomatrixPair::observations(uint j) const {
    return s.back().observations(j);
}

// The lowest min_genes setting
pair<size_t,size_t> PhenomatrixPair::min_genes() const {
    return make_pair<size_t,size_t>(p.top().min_genes(), s.back().min_genes());
}

id_set PhenomatrixPair::row_ids() const {
    return s.back().row_ids();
}

id_set PhenomatrixPair::predictable_column_ids() const {
    return p.top().column_ids();
}


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


// Create the predict matrix stack
stack<Phenomatrix> PhenomatrixPair::create_phenomatrix_stack(uint id, uint given_id, size_t min_genes) {
    stack<Phenomatrix> new_stack; new_stack.push(Phenomatrix(id, given_id, min_genes));
    return new_stack;
}

// Create the source matrix stack (well, actually a list)
list<PhenomatrixBase> PhenomatrixPair::create_phenomatrix_list(uint id, size_t min_genes) {
    list<PhenomatrixBase> new_list; new_list.push_back(PhenomatrixBase(id, true, min_genes));
    return new_list;
}
list<PhenomatrixBase> PhenomatrixPair::create_phenomatrix_list(const PhenomatrixBase& pb) {
    list<PhenomatrixBase> new_list; new_list.push_back(pb);
    return new_list;
}

#ifdef RICE
stack<Phenomatrix> PhenomatrixPair::create_phenomatrix_stack(Rice::Object predict_or_id, uint given_id, size_t min_genes) {
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

list<PhenomatrixBase> PhenomatrixPair::create_phenomatrix_list(Rice::Object source_or_id, size_t min_genes) {
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


// Cosine similarity between the two matrices in pair p, using TF-IDF, and
// treating both as a single corpus.
double cosine_similarity(const PhenomatrixPair* const p, uint j1, uint j2) {
    using boost::numeric::ublas::inner_prod;
    sparse_document_vectors v1v2 = p->document_vectors(j1, j2);

    double sim = inner_prod(v1v2.first, v1v2.second) / sqrt(inner_prod(v1v2.first, v1v2.first) * inner_prod(v1v2.second, v1v2.second));
    if (sim > 1.0) return 0.0;
    else return 1.0 - sim;
}

// Tanimoto coefficient between the two matrices in pair p, using TF-IDF, and
// treating both as a single corpus.
double tanimoto_coefficient(const PhenomatrixPair* const p, uint j1, uint j2) {
    using boost::numeric::ublas::inner_prod;
    sparse_document_vectors v1v2 = p->document_vectors(j1, j2);

    double sim = inner_prod(v1v2.first, v1v2.second) / ( inner_prod(v1v2.first, v1v2.first) + inner_prod(v1v2.second, v1v2.second) - inner_prod(v1v2.first, v1v2.second) );
    if (sim > 1.0) return 0.0;
    else return 1.0 - sim;
}

// Determine the manhattan magnitude of a vector (sum of its elements).
double manhattan_magnitude(const sparse_document_vector& v) {
    typedef sparse_document_vector::const_iterator sdv_iter;
    
    double accum = 0.0;
    for (sdv_iter it = v.begin(); it != v.end(); ++it)
        accum += std::abs(*it);

    return accum;
}

// Sample Pearson correlation coefficient.
double pearson(const PhenomatrixPair* const p, uint j1, uint j2) {
    using boost::numeric::ublas::inner_prod;
    sparse_document_vectors v1v2 = p->document_vectors(j1, j2);

    double xy = inner_prod(v1v2.first, v1v2.second);
    double xx = inner_prod(v1v2.first, v1v2.first);
    double yy = inner_prod(v1v2.second, v1v2.second);
    double magx = manhattan_magnitude(v1v2.first);
    double magy = manhattan_magnitude(v1v2.second);

    size_t n = p->max_intersection_size();

    return (n*xy - magx * magy) / sqrt((n*xx - magx*magx) * (n*yy - magy*magy));
}


double jaccard(const PhenomatrixPair* const p, uint j1, uint j2) {
    std::pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;

    if (p->min_idf() <= 0.0)
        // Simple: No idf threshold set, just go ahead and use the actual matrix.
        return jaccard(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );


    else {

        // Less simple: idf_threshold > 0.0. Compute document vectors and use those
        // as matrix instead.
        using boost::numeric::ublas::element_prod;

        tuple<size_t,size_t,size_t> mnk = sparse_drawn_defective_nnz( p->document_vectors(j1, j2) );
        if (mnk.get<0>() == 0 || mnk.get<1>() == 0 || mnk.get<2>() == 0) return 1.0;

        return jaccard(mnk.get<0>(), mnk.get<1>(), mnk.get<2>(), p->max_intersection_size() );
    }
}


double sorensen(const PhenomatrixPair* const p, uint j1, uint j2) {
    std::pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;

    if (p->min_idf() <= 0.0)
        // Simple: No idf threshold set, just go ahead and use the actual matrix.
        return sorensen(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );


    else {

        // Less simple: idf_threshold > 0.0. Compute document vectors and use those
        // as matrix instead.
        using boost::numeric::ublas::element_prod;

        tuple<size_t,size_t,size_t> mnk = sparse_drawn_defective_nnz( p->document_vectors(j1, j2) );
        if (mnk.get<0>() == 0 || mnk.get<1>() == 0 || mnk.get<2>() == 0) return 1.0;

        return sorensen(mnk.get<0>(), mnk.get<1>(), mnk.get<2>(), p->max_intersection_size() );
    }
}

// Exposes manhattan(m,n,k,N) and PhenomatrixPair to eachother.
double manhattan(const PhenomatrixPair* const p, uint j1, uint j2) {
    using boost::numeric::ublas::inner_prod;
    sparse_document_vectors v1v2 = p->document_vectors(j1, j2);
    sparse_document_vector v3(v1v2.first - v1v2.second);
    
    return manhattan_magnitude(v3) / (double)(p->max_intersection_size());
}

// Exposes euclidean(m,n,k,N) and PhenomatrixPair to eachother.
double euclidean(const PhenomatrixPair* const p, uint j1, uint j2) {
    using boost::numeric::ublas::inner_prod;
    sparse_document_vectors v1v2 = p->document_vectors(j1, j2);

    sparse_document_vector v3(v1v2.first - v1v2.second);
    return sqrt(inner_prod(v3, v3)) / (double)(p->max_intersection_size());
}

// Exposes hypergeometric(m,n,k,N) and PhenomatrixPair to eachother.
// Each distance function should have something like this -- and it's important
// that it have exactly these arguments.
double hypergeometric(const PhenomatrixPair* const p, uint j1, uint j2) {
    
    std::pair<size_t,size_t> obs_j1_j2 = p->observations_sizes(j1,j2);
    if (obs_j1_j2.first == 0 || obs_j1_j2.second == 0) return 1.0;
    
    if (p->min_idf() <= 0.0)
        // Simple: No idf threshold set, just go ahead and use the actual matrix.
        return hypergeometric(obs_j1_j2.first, obs_j1_j2.second, p->intersection(j1,j2).size(), p->max_intersection_size() );

    
    else {

        // Less simple: idf_threshold > 0.0. Compute document vectors and use those
        // as matrix instead.
        using boost::numeric::ublas::element_prod;
        
        tuple<size_t,size_t,size_t> mnk = sparse_drawn_defective_nnz( p->document_vectors(j1, j2) );
        if (mnk.get<0>() == 0 || mnk.get<1>() == 0 || mnk.get<2>() == 0) return 1.0;

        return hypergeometric(mnk.get<0>(), mnk.get<1>(), mnk.get<2>(), p->max_intersection_size() );
    }
}


size_t nnz(const sparse_document_vector& v) {
    size_t res = 0;
    for (sparse_document_vector::const_iterator it = v.begin(); it != v.end(); ++it)
        if (*it != 0.0) res++;

    return res;
}

tuple<size_t, size_t, size_t> sparse_drawn_defective_nnz(const sparse_document_vectors& v1v2) {
    size_t m = 0, n = 0, k = 0;
    // TODO: Speed this up by writing my own iterator, possibly?
    sparse_document_vector v_sum(v1v2.first + v1v2.second);
    for (sparse_document_vector::const_iterator it = v_sum.begin(); it != v_sum.end(); ++it) {
        if (v1v2.first[it.index()] != 0.0) {
            m++;
            if (v1v2.second[it.index()] != 0.0) {
                n++;
                k++;
            }
        } else {
            if (v1v2.second[it.index()] == 0.0) throw; // REMOVE THIS IF IT WORKS ONCE.
            n++;
        }
    }
    return make_tuple(m, n, k);
}