#include "distance_matrix.h"


// All functions which make use of Classifier must be located here since it's
// a circular include

DistanceMatrix::DistanceMatrix(conn_t* c_, uint predict_matrix_id, const id_set& source_matrix_ids, const string& distfn, const cparams& classifier_params)
 : c(c_),
   destroy_c(false),
   predict_matrix_(c, predict_matrix_id),
   distance_function(switch_distance_function(distfn))
{
    construct_classifier(classifier_params);

    for (id_set::const_iterator st = source_matrix_ids.begin(); st != source_matrix_ids.end(); ++st) {
        source_matrices.push_back( Phenomatrix(c, *st) );
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
    // Do not delete shared connections!
    if (destroy_c) delete c;
    delete classifier;
}


pcolumn DistanceMatrix::predict(uint j) const {
    return (*classifier)(j);
}


#ifdef RICE
using namespace Rice;

// This constructor is the one we use for the Ruby interface (RICE) since
// Ruby would likely have trouble with std::set. Instead, it takes an array.
DistanceMatrix::DistanceMatrix(const string& dbstr, uint predict_matrix_id, const Array& source_matrix_ids, const string& distfn, Rice::Object classifier_params_h)
 : c(new conn_t(dbstr)), destroy_c(true), predict_matrix_(c, predict_matrix_id),
   distance_function(switch_distance_function(distfn))
{
    // Convert from Hash to cparams
    cparams classifier_params( from_ruby<Rice::Symbol>(classifier_params_h.call("fetch", Rice::Symbol("classifier"))).str() );
    classifier_params.k = from_ruby<uint>(classifier_params_h.call("fetch", Rice::Symbol("k")));
    
    construct_classifier( classifier_params );

    for (Array::const_iterator st = source_matrix_ids.begin(); st != source_matrix_ids.end(); ++st) {
        uint id = from_ruby<uint>(*st);
#ifdef DEBUG_TRACE
        cerr << "distance_matrix.o: Adding phenomatrix " << id << " to distance matrix" << endl;
#endif
        source_matrices.push_back( Phenomatrix(c, id) );
    }
}



// IT IS CRITICAL THAT Rice:: TYPES LEAVE OFF THE NAMESPACE BEYOND THIS POINT!
// Remember, this is extern "C", and C doesn't understand namespaces!
// If you don't leave off the namespace, you'll get bizarre errors like "Can't convert from Hash to Rice::Hash."
extern "C"
void Init_distance_matrix() {

    Rice::Module rb_mFastknn = define_module("Fastknn");

    Data_Type<DistanceMatrix> rb_cDistanceMatrix =
            define_class_under<DistanceMatrix>(rb_mFastknn, "DistanceMatrix")
            .define_constructor(Constructor<DistanceMatrix,const string&, uint, Array, const string&, Object>())
            .define_method("source_matrix_ids", &DistanceMatrix::source_matrix_ids)
            .define_method("max_intersection_size", &DistanceMatrix::max_intersection_size)
            .define_method("max_intersection_count", &DistanceMatrix::max_intersection_size)
            .define_method("intersection_size", &DistanceMatrix::intersection_size)
            .define_method("intersection_count", &DistanceMatrix::intersection_size)
            .define_method("nearest_id", &DistanceMatrix::rb_nearest_id)
            .define_method("nearest_distance", &DistanceMatrix::rb_nearest_distance)
            .define_method("nearest", &DistanceMatrix::rb_nearest)
            .define_method("distance", &DistanceMatrix::rb_distance)
            .define_method("knearest", &DistanceMatrix::rb_knearest)
            .define_method("predict", &DistanceMatrix::rb_predict);
}

#endif

int main() {
    conn_t* c = new conn_t("dbname=crossval_development user=jwoods password=youwish1");
    id_set sources; sources.insert(3);
    cparams cp("naivebayes"); cp.k = 10;
    DistanceMatrix d(c, 185, sources, "hypergeometric", cp);
    cout << "Max common items: " << d.max_intersection_size() << endl;
    return 0;
}


