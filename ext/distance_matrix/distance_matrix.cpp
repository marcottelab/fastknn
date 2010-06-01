#include "distance_matrix.h"


// All functions which make use of Classifier must be located here since it's
// a circular include

DistanceMatrix::DistanceMatrix(
        uint predict_matrix_id,
        id_set source_matrix_ids,
        string distfn,
        cparams classifier_params
)
 : predict_matrix_(predict_matrix_id, predict_matrix_id)
{
    construct_classifier(classifier_params);

    for (id_set::const_iterator st = source_matrix_ids.begin(); st != source_matrix_ids.end(); ++st) {
        source_matrices.push_back( PhenomatrixPair(predict_matrix_id, *st, distfn) );
    }
}


void DistanceMatrix::construct_classifier(const cparams& classifier_params) {
    if (classifier_params.classifier == "naivebayes")
        classifier = new NaiveBayes(this, classifier_params.k);
    else {
        string err = "distance_matrix.o: Unrecognized classifier '" + classifier_params.classifier + "'!";
#ifdef RICE
        throw Rice::Exception(rb_eArgError, err.c_str());
#else
        cerr << err << endl;
        throw;
#endif
    }
}


DistanceMatrix::~DistanceMatrix() {
    delete classifier;
}


pcolumn DistanceMatrix::predict(uint j) const {
    return (*classifier)(j);
}


#ifdef RICE
using namespace Rice;

////////////////////////////////////////////////////////////////////////////////
// CONVERSION FUNCTIONS FOR RETURN VALUES AND ARGUMENT TYPES THAT NEED TO BE
// EXPOSED TO RUBY THROUGH RICE

template<>
Rice::Object to_ruby<id_dist>(id_dist const & d) {
    return to_ruby<Array>(d.to_a());
}

template<>
Rice::Object to_ruby<id_dist_iter>(id_dist_iter const & d) {
    return to_ruby<Array>(d.to_a());
}


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


template<>
Object to_ruby<cparams>(cparams const & param) {
    return param.to_h();
}

template<>
cparams from_ruby<cparams>(Object x) {
    Hash hash(x);
    cparams params( from_ruby<Symbol>(hash[ Symbol("classifier") ]).str() );
    params.k = from_ruby<uint>( hash[ Symbol("k") ]);

    return params;
}


template <>
Object to_ruby<proximity_queue>(proximity_queue const & d) {
    Array ary;
    proximity_queue pq(d); // make a copy
    while (pq.size() > 0) {
        ary.push(to_ruby<proximity_queue::value_type>(pq.top()));
        pq.pop();
    }
    return ary;
}


template <>
Object to_ruby<pcolumn>(pcolumn const & d) {
    Hash h;
    for (pcolumn::const_iterator i = d.begin(); i != d.end(); ++i) {
        h[to_ruby<uint>(i->first)] = to_ruby<float>(i->second);
    }
    return h;
}


// Convert from Rice::Array to std::set
template <>
pcolumn from_ruby<pcolumn>(Object x) {
    Hash h(x);
    pcolumn xmap;
    for (Hash::iterator i = h.begin(); i != h.end(); ++i)
        xmap[from_ruby<uint>(i->first)] = from_ruby<float>(i->second);
    return xmap;
}



// IT IS CRITICAL THAT Rice:: TYPES LEAVE OFF THE NAMESPACE BEYOND THIS POINT!
// Remember, this is extern "C", and C doesn't understand namespaces!
// If you don't leave off the namespace, you'll get bizarre errors like "Can't convert from Hash to Rice::Hash."
extern "C"
void Init_distance_matrix() {

    Rice::Module rb_mFastknn = define_module("Fastknn");

    #include "rice_connection.cpp"
    #include "rice_phenomatrix.cpp"

    Data_Type<DistanceMatrix> rb_cDistanceMatrix =
            define_class_under<DistanceMatrix>(rb_mFastknn, "DistanceMatrix")
            .define_constructor(Constructor<DistanceMatrix,uint,id_set,string,cparams>())
            .define_method("source_matrix_ids", &DistanceMatrix::source_matrix_ids)
            .define_method("intersection_size", &DistanceMatrix::intersection_size)
            .define_method("intersection_count", &DistanceMatrix::intersection_size)
            .define_method("nearest", &DistanceMatrix::nearest)
            .define_method("distance", &DistanceMatrix::distance)
            .define_method("knearest", &DistanceMatrix::knearest,
                           ( Arg("j"),
                             Arg("k") = (uint)(1),
                             Arg("bound") = (double)(1.0))       )
            .define_method("predict", &DistanceMatrix::predict)
            .define_method("predict_and_write", &DistanceMatrix::predict_and_write,
                           ( Arg("j"),
                             Arg("write_rows") = id_set())     )
            .define_method("predict_and_write_all", &DistanceMatrix::predict_and_write_all,
                           ( Arg("write_rows") = id_set())     )
            ;
}


#else

int main() {
    id_set sources; sources.insert(3);
    cparams cp("naivebayes"); cp.k = 10;
    DistanceMatrix d(185, sources, "hypergeometric", cp);
    return 0;
}

#endif
