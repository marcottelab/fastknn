    Data_Type<Connection> rb_cConnection =
            define_class_under<Connection>(rb_mFastknn, "Connection")
            .define_constructor(Constructor<Connection>())
            .define_method("connect", &Connection::connect)
            .define_method("connected?", &Connection::connected)
            .define_method("test_singleton", &Connection::instance)
            .define_method("count", &Connection::count)
            .define_method("instance", &Connection::instance);