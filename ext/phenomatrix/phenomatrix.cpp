#include "phenomatrix_base.h"
#include "phenomatrix.h"


#ifdef RICE
using namespace Rice;

extern "C"
void Init_phenomatrix() {

    Rice::Module rb_mFastknn = define_module("Fastknn");

    Data_Type<PhenomatrixBase> rb_cPhenomatrixBase =
            define_class_under<PhenomatrixBase>(rb_mFastknn, "PhenomatrixBase")
            .define_constructor(Constructor<PhenomatrixBase,const std::string&, uint>())
            .define_method("parent_id", &PhenomatrixBase::rb_parent_id)
            .define_method("root_id", &PhenomatrixBase::rb_root_id)
            .define_method("id", &PhenomatrixBase::id)
            .define_method("observations_count", &PhenomatrixBase::observations_size)
            .define_method("has_column?", &PhenomatrixBase::has_column)
            .define_method("row_count", &PhenomatrixBase::row_count);

    Data_Type<Phenomatrix> rb_cPhenomatrix =
            define_class_under<Phenomatrix, PhenomatrixBase>(rb_mFastknn, "Phenomatrix")
            .define_constructor(Constructor<Phenomatrix,const std::string&, uint, uint>())
            .define_method("source_id", &Phenomatrix::source_id);
}
#endif

//#ifndef DISTANCE_MATRIX_H_
//int main() {
//    Phenomatrix p("dbname=crossval_development user=jwoods password=youwish1", 185);
//    return 0;
//}
//#endif

