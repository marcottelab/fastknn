#ifndef CPARAMS_H_
# define CPARAMS_H_

#include <string>

#ifdef RICE
#include <rice/Hash.hpp>
#include <rice/Symbol.hpp>
#endif

using std::string;

class cparams {
public:
    cparams(const string& name) : classifier(name), k(0), max_distance(1.0) { }

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