#include "connection.h"

string join(const id_set& ids, const string& join_str) {
    list<string> str_ids;
    BOOST_FOREACH( uint id, ids )
    {
        str_ids.push_back(lexical_cast<string>(id));
    }
    return boost::algorithm::join(str_ids, join_str);
}

// Fetches an unsigned integer type (like a count) which will never be nil.
uint Connection::fetch_count(const string& sql) {
    work_t w(*c);
    result_t r = w.exec(sql);

    uint ret;
    r[0][0].to(ret);
    return ret;
}


string Connection::fetch_type(const string& table, uint id) {
    work_t w(*c);
    result_t r = w.exec("SELECT type FROM " + table + " WHERE id = " + lexical_cast<string>(id) + ";");

    return string(r[0][0].c_str());
}

// Note: Make sure to include ORDER BY whatever the ID column is. Let SQL do the
// sorting.
id_set Connection::fetch_id_set(const string& sql) {
    work_t w(*c);
    result_t r = w.exec(sql);

    id_set ret;
    id_set::iterator hint = ret.end();

    for (result_t::const_iterator rt = r.begin(); rt != r.end(); ++rt) {
        uint id = 0;
        (*rt)[0].to(id);
        if (id != 0) hint = ret.insert(hint, id);
    }

    return ret;
}


// Just like fetch_count, but checks for nil
size_t Connection::fetch_id(const string& sql) {
    work_t w(*c);
    result_t r = w.exec(sql);

    if (r[0][0].is_null())
        return 0; // 0 indicates NULL since no ID will ever be 0
    else {
        uint ret;
        r[0][0].to(ret);
        return ret;
    }
}

conn_t* Connection::c = NULL;
size_t  Connection::count_ = 0;

#ifdef RICE
#include <rice/Data_Type.hpp>
#include <rice/Constructor.hpp>
#include <rice/Module.hpp>

using namespace Rice;

extern "C"
void Init_connection() {

    Module rb_mFastknn = define_module("Fastknn");

    #include "rice_connection.cpp"
}

#endif

