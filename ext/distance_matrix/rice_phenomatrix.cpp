    Data_Type<PhenomatrixBase> rb_cPhenomatrixBase =
            define_class_under<PhenomatrixBase>(rb_mFastknn, "PhenomatrixBase")
            .define_constructor(Constructor<PhenomatrixBase,uint>())
            .define_method("parent_id", &PhenomatrixBase::rb_parent_id)
            .define_method("root_id", &PhenomatrixBase::rb_root_id)
            .define_method("id", &PhenomatrixBase::id)
            .define_method("observations", &PhenomatrixBase::observations)
            .define_method("observations_count", &PhenomatrixBase::observations_size)
            .define_method("has_column?", &PhenomatrixBase::has_column)
            .define_method("row_count", &PhenomatrixBase::row_count)
            .define_method("row_ids", &PhenomatrixBase::row_ids)
            .define_method("column_ids", &PhenomatrixBase::column_ids)
            .define_method("child_ids", &PhenomatrixBase::child_ids)
            .define_method("child_row_ids", &PhenomatrixBase::child_row_ids)
            ;

    Data_Type<Phenomatrix> rb_cPhenomatrix =
            define_class_under<Phenomatrix, PhenomatrixBase>(rb_mFastknn, "Phenomatrix")
            .define_constructor(Constructor<Phenomatrix, uint, uint>())
            .define_method("source_id", &Phenomatrix::source_id);

    Data_Type<FusionPhenomatrix> rb_cFusionPhenomatrix =
            define_class_under<FusionPhenomatrix, PhenomatrixBase>(rb_mFastknn, "FusionPhenomatrix")
            .define_constructor(Constructor<FusionPhenomatrix, uint, id_set>())
            .define_method("source_ids", &FusionPhenomatrix::source_ids);

    // Note that the constructor has an extra object to accomodate the self bug.
    // This is never used, from my perspective.
    Data_Type<PhenomatrixPair> rb_cPhenomatrixPair =
            define_class_under<PhenomatrixPair>(rb_mFastknn, "PhenomatrixPair")
            .define_constructor(Constructor<PhenomatrixPair, Object, Object, Object>())
            .define_method("distance", &PhenomatrixPair::distance)
            .define_method("nearest", &PhenomatrixPair::nearest)
            .define_method("predict_matrix_has_column?", &PhenomatrixPair::predict_matrix_has_column)
            .define_method("source_matrix_has_column?", &PhenomatrixPair::source_matrix_has_column)
            .define_method("push_mask", &PhenomatrixPair::push_mask)
            .define_method("pop_mask", &PhenomatrixPair::pop_mask)
            .define_method("distance_function=", &PhenomatrixPair::set_distance_function)
            ;