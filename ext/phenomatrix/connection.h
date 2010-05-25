#ifndef CONNECTION_H_
# define CONNECTION_H_

#include <pqxx/connection.hxx>
#include <pqxx/transaction.hxx>
#include <pqxx/result.hxx>

#include <string>
#include <set>
#include <boost/lexical_cast.hpp>

using boost::lexical_cast;
using std::string;
typedef pqxx::result result_t;
typedef pqxx::connection conn_t;
typedef pqxx::work work_t;
typedef std::set<uint> id_set;

uint fetch_count(conn_t* c, const string& sql);
string fetch_type(conn_t* c, const string& table, uint id);
id_set fetch_id_set(conn_t* c, const string& sql);
size_t fetch_id(conn_t* c, const string& sql);

#endif
