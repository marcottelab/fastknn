#include "connection.h"

// Fetches an unsigned integer type (like a count) which will never be nil.
uint fetch_count(conn_t* c, const string& sql) {
    work_t w(*c);
    result_t r = w.exec(sql);

    uint ret;
    r[0][0].to(ret);
    return ret;
}


string fetch_type(conn_t* c, const string& table, uint id) {
    work_t w(*c);
    result_t r = w.exec("SELECT type FROM " + table + " WHERE id = " + lexical_cast<string>(id) + ";");

    return string(r[0][0].c_str());
}

// Note: Make sure to include ORDER BY whatever the ID column is. Let SQL do the
// sorting.
id_set fetch_id_set(conn_t* c, const string& sql) {
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
size_t fetch_id(conn_t* c, const string& sql) {
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