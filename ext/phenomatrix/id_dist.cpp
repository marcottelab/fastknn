#include "id_dist.h"

#ifdef RICE
Rice::Array id_dist_iter::to_a() const {
    Rice::Array a;
    a.push(to_ruby<uint>(id));
    a.push(to_ruby<double>(distance));

    // push the source matrix ID
    if (matrix_id > 0)
        a.push(to_ruby<uint>(matrix_id));
    else
        a.push(to_ruby<uint>(matrix_iter->id()));

    return a;
}


template <typename T, typename U>
Rice::Array pair_to_array(const std::pair<T,U>& vals) {
    Rice::Array a;
    a.push( to_ruby<T>(vals.first) );
    a.push( to_ruby<U>(vals.second) );
    return a;
}

#endif

id_dist_iter::id_dist_iter(const uint& id_, const double& d_, matrix_list::const_iterator m_iter_, uint matrix_id_)
    : id_dist(id_, d_), matrix_iter(m_iter_), matrix_id(matrix_id_)
{ }

id_dist_iter::id_dist_iter(const uint& id_, const double& d_, uint matrix_id_)
    : id_dist(id_, d_), matrix_id(matrix_id_)
{ }

