    Data_Type<DistanceMatrix> rb_cDistanceMatrix =
            define_class_under<DistanceMatrix>(rb_mFastknn, "DistanceMatrix")
            .define_constructor(Constructor<DistanceMatrix, uint, Object, size_t>())
            .define_method("id", &DistanceMatrix::id)
            .define_method("source_matrix_ids", &DistanceMatrix::source_matrix_ids)
            .define_method("intersection_size", &DistanceMatrix::intersection_size)
            .define_method("intersection_count", &DistanceMatrix::intersection_size)
            .define_method("nearest", &DistanceMatrix::nearest)
            .define_method("distance", &DistanceMatrix::distance)
            .define_method("distance_function=", &DistanceMatrix::set_distance_function)
            .define_method("distance_functions", &DistanceMatrix::get_distance_functions)
            .define_method("min_idf=", &DistanceMatrix::set_min_idf)
            .define_method("min_idfs", &DistanceMatrix::min_idfs)
            .define_method("knearest", &DistanceMatrix::knearest,
                           ( Arg("j"),
                             Arg("k") = (uint)(1),
                             Arg("bound") = (double)(1.0))     )
            .define_method("classifier=", &DistanceMatrix::set_classifier)
            .define_method("classifier", &DistanceMatrix::get_classifier)
            .define_method("predict", &DistanceMatrix::predict)
            .define_method("predict_and_write", &DistanceMatrix::predict_and_write,
                           ( Arg("j"),
                             Arg("write_rows") = id_set())     )
            .define_method("predict_and_write_all", &DistanceMatrix::predict_and_write_all,
                           ( Arg("write_rows") = id_set())     )
            .define_method("predict_and_write_to", &DistanceMatrix::predict_and_write_to,
                           ( Arg("dir"),
                             Arg("j"),
                             Arg("write_rows"))                )
            .define_method("push_mask", &DistanceMatrix::push_mask)
            .define_method("pop_mask", &DistanceMatrix::pop_mask)
            .define_method("predict_matrix_has_column?", &DistanceMatrix::predict_matrix_has_column)
            .define_method("predictable_columns", &DistanceMatrix::predictable_columns)
            .define_method("predict_matrix", &DistanceMatrix::predict_matrix)
            .define_method("source_matrix_pairs", &DistanceMatrix::source_matrix_pairs)
            .define_method("crossvalidate", &DistanceMatrix::crossvalidate)
            .define_method("min_genes", &DistanceMatrix::min_genes)
            ;
