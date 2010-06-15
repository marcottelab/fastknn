#ifndef PARAMS_H_
# define PARAMS_H_

#include <string>

#ifdef RICE
#include <rice/Hash.hpp>
#include <rice/Symbol.hpp>
#endif

using std::string;


class params {
public:
    params() { }
    virtual ~params() { }
    
    virtual Rice::Object to_h() const = 0;
};


class classifier_params : params {
public:
    classifier_params(const string& name) : classifier(name), k(0), max_distance(1.0) { }
    classifier_params() : classifier(), k(0), max_distance(1.0) { }

    bool operator==(const classifier_params& rhs) const {
        return (classifier == rhs.classifier && k == rhs.k && max_distance == rhs.max_distance);
    }

    bool operator!=(const classifier_params& rhs) const {
        return (classifier != rhs.classifier || k != rhs.k || max_distance != rhs.max_distance);
    }

    string classifier;
    uint k;
    float max_distance;

#ifdef RICE
    Rice::Object to_h() const {
        Rice::Hash h;

        h[ Rice::Symbol("classifier")   ]   = to_ruby<string>(classifier);
        h[ Rice::Symbol("k")            ]   = to_ruby<uint>(k);
        h[ Rice::Symbol("max_distance") ]   = to_ruby<float>(max_distance);
        
        return h;
    }
#endif
};


#endif