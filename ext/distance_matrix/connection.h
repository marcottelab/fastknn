#ifndef CONNECTION_H_
# define CONNECTION_H_

#ifdef RICE
#include <rice/Exception.hpp>
#endif

#include <boost/shared_ptr.hpp>

#include <pqxx/connection.hxx>
#include <pqxx/transaction.hxx>
#include <pqxx/result.hxx>

#include <map>
#include <utility>
#include <string>
#include <set>
#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/foreach.hpp>
#include <list>


using std::list;
using std::map;
using std::make_pair;
using std::cerr;
using std::endl;
using boost::shared_ptr;
using boost::lexical_cast;
using std::string;
typedef pqxx::result result_t;
typedef pqxx::connection conn_t;
typedef pqxx::work work_t;
typedef std::set<uint> id_set;


string join(const id_set& ids, const string& join_str);

class Connection {
public:
    Connection() {
        ++count_;
    }

    ~Connection() {
        --count_;
        if (count_ == 0 && c != NULL) {
            delete c;
            c = NULL;
        }
    }

    size_t count() const { return count_; }

    static Connection& instance() {
        static Connection the_connection;
#ifdef DEBUG_CONNECTION_POINTER
        cerr << "c: '" << c << "'" << endl;
#endif
        return the_connection;
    }

    bool connected() const { return c != 0; }

    // Attempt to connect unless we've already connected.
    // Return true if there's a connection after the effort completes, otherwise
    // return false.
    bool connect(const string& dbstr) {
        if (c == NULL) c = new conn_t(dbstr);
        if (c)  return true;
        return false;
    }


    work_t* work() {
        if (connected())
            return new work_t(*c);
        else
#ifdef RICE
            throw Rice::Exception(rb_eArgError, "connection.h: work(): Not connected to database!");
#else
            throw;
#endif
    }


    uint fetch_count(const string& sql);
    string fetch_type(const string& table, uint id);
    id_set fetch_id_set(const string& sql);
    size_t fetch_id(const string& sql);

    template <typename K, typename V>
    map<K,V> fetch_map(const string& sql) {
        work_t w(*c);
        result_t r = w.exec(sql);

        map<K,V> ret;

        typename map<K,V>::iterator hint(ret.end());

        for (result_t::const_iterator rt = r.begin(); rt != r.end(); ++rt) {
            K id = 0;
            V t;
            (*rt)[0].to(id);
            (*rt)[1].to(t);
            if (id != 0)
                hint = ret.insert(hint,
                                  make_pair( id, t )
                                 );
        }

        return ret;
    }

private:    
    static conn_t* c;
    static size_t count_;
};


#endif // CONNECTION_H_
