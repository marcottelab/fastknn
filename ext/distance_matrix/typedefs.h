#ifndef TYPEDEFS_H_
# define TYPEDEFS_H_

#ifdef RICE
#include <rice/Object.hpp>
#include <rice/Array.hpp>
#include <rice/Hash.hpp>
#include <rice/Module.hpp>
#include <rice/Data_Type.hpp>
#include <rice/Constructor.hpp>
using Rice::Array;
using Rice::Object;
using Rice::Hash;
#endif

#include <boost/unordered_map.hpp>
#include <utility>
#include <list>
#include <vector>
#include <queue>

#include "phenomatrix.h"

typedef std::list<Phenomatrix>              matrix_list;
typedef std::pair<uint,double>              id_dist_pair;
typedef std::list<id_dist_pair>             proximity_list;
// typedef type_shield<double,pair<uint, matrix_list::const_iterator> >            dist_id;
typedef boost::unordered_map<uint, float>   pcolumn; // prediction column

// This is the return value for
class id_dist_iter {
public:
    id_dist_iter(const uint& id_, const double& d_, matrix_list::const_iterator m_iter_)
    : id(id_), distance(d_), matrix_iter(m_iter_)
    { }

    uint  id;
    double distance;
    matrix_list::const_iterator matrix_iter;

    bool operator<(const id_dist_iter& rhs) const {
        return distance < rhs.distance;
    }

    bool operator>(const id_dist_iter& rhs) const {
        return distance > rhs.distance;
    }

#ifdef RICE
    Rice::Array to_a() const {
        Rice::Array a;
        a.push(to_ruby<uint>(id));
        a.push(to_ruby<double>(distance));
        a.push(to_ruby<uint>(matrix_iter->id())); // push the source matrix ID
        return a;
    }
#endif
};


#ifdef RICE
template <typename T, typename U>
Rice::Array pair_to_array(const std::pair<T,U>& vals) {
    Rice::Array a;
    a.push( to_ruby<T>(vals.first) );
    a.push( to_ruby<U>(vals.second) );
    return a;
}
#endif

#endif //TYPEDEFS_H_