#include "phenomatrix.h"

#ifdef RICE
using namespace Rice;

extern "C"
void Init_phenomatrix() {

    Rice::Module rb_mFastknn = define_module("Fastknn");

    Data_Type<Phenomatrix> rb_cPhenomatrix =
            define_class_under<Phenomatrix>(rb_mFastknn, "Phenomatrix")
            .define_constructor(Constructor<Phenomatrix,const std::string&, uint>())
            .define_method("parent_id", &Phenomatrix::rb_parent_id)
            .define_method("root_id", &Phenomatrix::rb_root_id)
            .define_method("id", &Phenomatrix::id)
            .define_method("observations_count", &Phenomatrix::observations_size)
            .define_method("has_column?", &Phenomatrix::has_column)
            .define_method("naive_row_count", &Phenomatrix::naive_row_count);
}
#endif

//#ifndef DISTANCE_MATRIX_H_
//int main() {
//    Phenomatrix p("dbname=crossval_development user=jwoods password=youwish1", 185);
//    return 0;
//}
//#endif

