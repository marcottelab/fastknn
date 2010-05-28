#include "phenomatrix_base.h"
#include "phenomatrix.h"
#include "phenomatrix_pair.h"


// Find the single nearest neighbor (FIRST FOUND, not equivalence class)
id_dist_iter PhenomatrixPair::nearest(uint j) const {

    double min_dist = 100;
    uint min_dist_id = 0;

    id_set s_columns = s.top().column_ids();

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
    const id_set& columns = s.top().column_ids();

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


#ifdef RICE

using namespace Rice;

template <>
Object to_ruby<id_set>(id_set const & d) {
    Array ary;
    for (id_set::const_iterator i = d.begin(); i != d.end(); ++i)
        ary.push( to_ruby<uint>(*i) );
    return ary;
}



// Convert from Rice::Array to std::set
template <>
id_set from_ruby<id_set >(Object x) {
    Array ary(x);
    id_set result;
    for (Array::iterator i = ary.begin(); i != ary.end(); ++i)
        result.insert(from_ruby<uint>(*i));
    return result;
}


extern "C"
void Init_phenomatrix() {

    Rice::Module rb_mFastknn = define_module("Fastknn");

    Data_Type<Connection> rb_cConnection =
            define_class_under<Connection>(rb_mFastknn, "Connection")
            .define_constructor(Constructor<Connection>())
            .define_method("connect", &Connection::connect)
            .define_method("connected?", &Connection::connected)
            .define_method("test_singleton", &Connection::instance)
            .define_method("count", &Connection::count)
            .define_method("instance", &Connection::instance);

    Data_Type<PhenomatrixBase> rb_cPhenomatrixBase =
            define_class_under<PhenomatrixBase>(rb_mFastknn, "PhenomatrixBase")
            .define_constructor(Constructor<PhenomatrixBase,uint>())
            .define_method("parent_id", &PhenomatrixBase::rb_parent_id)
            .define_method("root_id", &PhenomatrixBase::rb_root_id)
            .define_method("id", &PhenomatrixBase::id)
            .define_method("observations_count", &PhenomatrixBase::observations_size)
            .define_method("has_column?", &PhenomatrixBase::has_column)
            .define_method("row_count", &PhenomatrixBase::row_count);

    Data_Type<Phenomatrix> rb_cPhenomatrix =
            define_class_under<Phenomatrix, PhenomatrixBase>(rb_mFastknn, "Phenomatrix")
            .define_constructor(Constructor<Phenomatrix, uint, uint>())
            .define_method("source_id", &Phenomatrix::source_id);

    Data_Type<PhenomatrixPair> rb_cPhenomatrixPair =
            define_class_under<PhenomatrixPair>(rb_mFastknn, "PhenomatrixPair")
            .define_constructor(Constructor<PhenomatrixPair, uint, uint, const string&>())
            .define_method("distance", &PhenomatrixPair::distance)
            .define_method("nearest", &PhenomatrixPair::nearest)
            .define_method("predict_matrix_has_column?", &PhenomatrixPair::predict_matrix_has_column)
            .define_method("source_matrix_has_column?", &PhenomatrixPair::source_matrix_has_column)
            .define_method("push_mask", &PhenomatrixPair::push_mask)
            .define_method("pop_mask", &PhenomatrixPair::pop_mask)
            ;

}
#else
int main() {
    PhenomatrixPair p(185, 3, "hypergeometric");
    return 0;
}
#endif

//#ifndef DISTANCE_MATRIX_H_
//int main() {
//    Phenomatrix p("dbname=crossval_development user=jwoods password=youwish1", 185);
//    return 0;
//}
//#endif

