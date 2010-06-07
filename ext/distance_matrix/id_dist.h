#ifndef ID_DIST_H_
# define ID_DIST_H_

#include "typedefs.h"
#include "phenomatrix_pair.h"

class PhenomatrixPair;
typedef std::list<PhenomatrixPair>          matrix_list;

class id_dist {
public:
    id_dist(const uint& id_, const double& d_)
    : id(id_), distance(d_)
    { }

    bool operator<(const id_dist& rhs) const {
        return distance < rhs.distance;
    }

    bool operator>(const id_dist& rhs) const {
        return distance > rhs.distance;
    }

    uint  id;
    double distance;

#ifdef RICE
    virtual Rice::Array to_a() const {
        Rice::Array a;
        a.push(to_ruby<uint>(id));
        a.push(to_ruby<double>(distance));
        return a;
    }
#endif
};

// This is the return value for
class id_dist_iter : public id_dist {
public:
    id_dist_iter(const uint& id_, const double& d_, matrix_list::const_iterator m_iter_, uint matrix_id_ = 0);
    id_dist_iter(const uint& id_, const double& d_, uint matrix_id_);
    matrix_list::const_iterator matrix_iter;
    uint matrix_id;

#ifdef RICE
    Rice::Array to_a() const;
#endif
};


#endif