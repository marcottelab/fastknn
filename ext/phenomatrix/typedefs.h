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

const double MAX_DISTANCE = 1.0;

typedef std::pair<uint,double>              id_dist_pair;
typedef std::list<id_dist_pair>             proximity_list;
typedef boost::unordered_map<uint, float>   pcolumn; // prediction column

#endif //TYPEDEFS_H_