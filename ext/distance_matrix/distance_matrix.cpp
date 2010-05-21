#include "distance_matrix.h"

#ifdef RICE
using namespace Rice;

extern "C"
void Init_distance_matrix() {

    Data_Type<DistanceMatrix> rb_cDistanceMatrix =
            define_class<DistanceMatrix>("DistanceMatrix")
            .define_constructor(Constructor<DistanceMatrix,const string&, uint, uint, string>())
            .define_method("max_intersection_size", &DistanceMatrix::max_intersection_size)
            .define_method("max_intersection_count", &DistanceMatrix::max_intersection_size)
            .define_method("intersection_size", &DistanceMatrix::intersection_size)
            .define_method("intersection_count", &DistanceMatrix::intersection_size)
            .define_method("nearest", &DistanceMatrix::rb_nearest)
            .define_method("nearest_distance", &DistanceMatrix::rb_nearest_distance)
            .define_method("distance", &DistanceMatrix::distance);
}

#endif

int main() {
    DistanceMatrix d("dbname=crossval_development user=jwoods password=youwish1", 185, 193, "hypergeometric");
    cout << "Max common items: " << d.max_intersection_size() << endl;
    return 0;
}



