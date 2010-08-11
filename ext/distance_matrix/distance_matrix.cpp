#include "distance_matrix.h"


// All functions which make use of Classifier must be located here since it's
// a circular include

#ifdef RICE
#include "ruby.h"

PhenomatrixPair DistanceMatrix::construct_source_matrix_pair(uint predict_matrix_id, Rice::Object source_or_id, size_t min_genes) {
    using namespace Rice;

    if (source_or_id.is_a( rb_cFixnum ))
        return PhenomatrixPair(predict_matrix_id, from_ruby<uint>(source_or_id), min_genes);

    if (source_or_id.is_a(Data_Type<PhenomatrixPair>::klass()))
        return from_ruby<PhenomatrixPair>(source_or_id);

    if (source_or_id.is_a(Data_Type<PhenomatrixBase>::klass()))
        return PhenomatrixPair(predict_matrix_id, from_ruby<PhenomatrixBase>(source_or_id), min_genes);

    throw Rice::Exception(rb_eArgError, "distance_matrix.cpp: construct_source_matrix_pair: Need a positive Fixnum, a PhenomatrixPair, or a PhenomatrixBase");
}

matrix_list DistanceMatrix::construct_source_matrices(uint predict_matrix_id, Rice::Object sources_or_ids, size_t min_genes) {
    using namespace Rice;

    matrix_list sources_list;
    if (sources_or_ids.is_a( rb_cArray )) {

        Rice::Array ary(sources_or_ids);
        for (Rice::Array::iterator it = ary.begin(); it != ary.end(); ++it)
            sources_list.push_back( construct_source_matrix_pair(predict_matrix_id, *it, min_genes) );

    }
    else // Not a list -- just one
        sources_list.push_back(construct_source_matrix_pair(predict_matrix_id, sources_or_ids, min_genes));
        
    return sources_list;
}

DistanceMatrix::DistanceMatrix(
        uint predict_matrix_id,
        Rice::Object sources_or_ids,
        size_t min_genes
)
 : source_matrices(construct_source_matrices(predict_matrix_id, sources_or_ids, min_genes)),
   predict_matrix_(predict_matrix_id, source_matrices),
   classifier_parameters("naivebayes"),
   classifier(NULL),
   min_genes_(min_genes)
{
    construct_classifier(classifier_parameters);
}


#else
DistanceMatrix::DistanceMatrix(
        uint predict_matrix_id,
        id_set source_matrix_ids,
        size_t min_genes
)
 : source_matrices(construct_source_matrices(predict_matrix_id, source_matrix_ids, min_genes)),
   predict_matrix_(predict_matrix_id, source_matrices),
   classifier_parameters("naivebayes"),
   classifier(NULL),
   min_genes_(min_genes)
{
    construct_classifier(classifier_parameters);
}
#endif

DistanceMatrix::DistanceMatrix(const DistanceMatrix& rhs)
: source_matrices(rhs.source_matrices),
  predict_matrix_(rhs.predict_matrix_),
  classifier_parameters(rhs.classifier_parameters),
  min_genes_(rhs.min_genes_)
{
    // Only initialize the classifier if the RHS DistanceMatrix has one setup.
    if (rhs.classifier)     construct_classifier(classifier_parameters);
    else                    classifier = NULL;
}


void DistanceMatrix::construct_classifier(const classifier_params& classifier_params) {
    if (classifier) delete classifier;
    
    classifier_parameters = classifier_params;

    if (classifier_params.classifier == "naivebayes")
        classifier = new NaiveBayes(this, classifier_params.k, classifier_params.max_distance, classifier_params.distance_exponent);
    else if (classifier_params.classifier == "average")
        classifier = new AverageClassifier(this, classifier_params.k, classifier_params.max_distance);
    else if (classifier_params.classifier == "simple")
        classifier = new SimpleClassifier(this);
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
    throw_on_missing_predict_column(j);
    
    return (*classifier)(j);
}


#ifdef RICE

#include "ruby_conversions.cpp"

Rice::Object DistanceMatrix::get_classifier() const {
    return to_ruby<classifier_params>(classifier_parameters);
}

using namespace Rice;




// IT IS CRITICAL THAT Rice:: TYPES LEAVE OFF THE NAMESPACE BEYOND THIS POINT!
// Remember, this is extern "C", and C doesn't understand namespaces!
// If you don't leave off the namespace, you'll get bizarre errors like "Can't convert from Hash to Rice::Hash."
extern "C"
void Init_distance_matrix() {

    Rice::Module rb_mFastknn = define_module("Fastknn");

    #include "rice_connection.cpp"
    #include "rice_phenomatrix.cpp"
    #include "rice_distance_matrix.cpp"
}


#else

int main() {
    Connection c;
    c.connect("dbname=crossval_development user=jwoods password=youwish1");

    id_set sources; sources.insert(3);
    classifier_params cp("naivebayes"); cp.k = 10;
    DistanceMatrix d(185, sources, "hypergeometric", cp);
    d.crossvalidate();
    return 0;
}

#endif
