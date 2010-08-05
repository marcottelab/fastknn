using namespace Rice;

////////////////////////////////////////////////////////////////////////////////
// CONVERSION FUNCTIONS FOR RETURN VALUES AND ARGUMENT TYPES THAT NEED TO BE
// EXPOSED TO RUBY THROUGH RICE

Data_Type<PhenomatrixPair> phenomatrix_pair_type;
//namespace Rice {
//    template <>
//    struct Default_Allocation_Strategy<PhenomatrixBase> {
//        static void free(PhenomatrixBase * obj) { }
//    };
//
////    template <>
////    struct Default_Allocation_Strategy<PhenomatrixPair> {
////        static void free(PhenomatrixPair * obj) { }
////    };
//}
//template<>
//PhenomatrixPair from_ruby<PhenomatrixPair>(Rice::Object x) {
//    Rice::Data_Object<PhenomatrixPair> d(x, phenomatrix_pair_type);
//    return *d;
//}

template<>
Rice::Object to_ruby<pair<size_t, size_t> >(pair<size_t, size_t> const & d) {
    Rice::Array ary;
    ary.push(d.first);
    ary.push(d.second);
    return ary;
}

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
id_set from_ruby<id_set>(Object x) {
    Array ary(x);
    id_set result;
    for (Array::iterator i = ary.begin(); i != ary.end(); ++i)
        result.insert(from_ruby<uint>(*i));
    return result;
}

template <>
Rice::Object to_ruby<fs::path>(fs::path const & d) {
    String s(d.string());
    return s;
}

template<>
Rice::Object to_ruby<map<uint, id_set> >(map<uint,id_set> const & d) {
    Hash h;
    for (map<uint,id_set>::const_iterator dt = d.begin(); dt != d.end(); ++dt)
        h[ to_ruby<uint>(dt->first) ] = to_ruby<id_set>(dt->second);
    return h;
}

template <>
Rice::Object to_ruby<set<fs::path> >(set<fs::path> const & d) {
    Array ary;
    for (set<fs::path>::const_iterator dt = d.begin(); dt != d.end(); ++dt)
        ary.push( to_ruby<fs::path>(*dt) );
    return ary;
}


template<>
Object to_ruby<classifier_params>(classifier_params const & param) {
    return param.to_h();
}

// All parameters must be included in the hash! If any are left out, this will
// probably throw an exception.
template<>
classifier_params from_ruby<classifier_params>(Object x) {
    Hash hash(x);

    const Symbol CLASSIFIER_S("classifier");
    const Symbol K_S("k");
    const Symbol MAX_DISTANCE_S("max_distance");
    const Symbol DISTANCE_EXPONENT_S("distance_exponent");

    classifier_params p(  from_ruby<Symbol>(hash[ CLASSIFIER_S ]).str() );
    p.k                 = from_ruby<uint>( hash[ K_S ]);
    p.max_distance      = from_ruby<float>( hash[ MAX_DISTANCE_S ] );
    p.distance_exponent = from_ruby<float>( hash[ DISTANCE_EXPONENT_S ] );

    return p;
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


template <>
Object to_ruby<sparse_document_vector>(sparse_document_vector const & d) {
    Hash h;
    for (sparse_document_vector::const_iterator i = d.begin(); i != d.end(); ++i) {
        h[to_ruby<uint>(i.index())] = to_ruby<double>(*i);
    }
    return h;
}