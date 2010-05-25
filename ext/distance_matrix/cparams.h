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
    cparams(const string& name) : classifier(name), k(0) { }

    string classifier;
    uint k;

#ifdef RICE
    Rice::Object to_h() const {
        Rice::Hash h;

        h[ Rice::Symbol("classifier")   ]   = to_ruby<string>(classifier);
        if (k > 0) h[ Rice::Symbol("k") ]   = to_ruby<uint>(k);
        
        return h;
    }
#endif
};


#endif