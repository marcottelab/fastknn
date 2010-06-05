#ifndef FUSION_PHENOMATRIX_H_
# define FUSION_PHENOMATRIX_H_

#include "phenomatrix_base.h"

#include <sstream>
#include <string>

using std::string;
using std::ostringstream;


class FusionPhenomatrix : public PhenomatrixBase {
public:
    FusionPhenomatrix(uint id, id_set given_ids)
    : PhenomatrixBase(id, true),
            given_ids_(given_ids)
    {
        if (given_ids.size() > 1 || *(given_ids.begin()) != id)
            inherit_construct();
    }

    FusionPhenomatrix(const FusionPhenomatrix& rhs)
    : PhenomatrixBase(rhs), given_ids_(rhs.given_ids_)
    { }

    FusionPhenomatrix(const FusionPhenomatrix& rhs, const id_set& mask_rows)
    : PhenomatrixBase(rhs, mask_rows), given_ids_(rhs.given_ids_)
    { }

    ~FusionPhenomatrix() { }

    // Get source matrix IDs
    id_set source_ids() const { return given_ids_; }

protected:


    // OVERRIDE SQL FROM PHENOMATRIX_BASE.
    string load_matrix_sql(uint matrix_id) const {
        ostringstream sql;
        sql   << "SELECT DISTINCT e.i,e.j FROM entries e " << endl
              << "INNER JOIN entries es ON (e.i = es.i) " << endl
              << "WHERE e.matrix_id = " << matrix_id << endl
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ")" << endl
              << " AND  e.type = 'Cell';" << endl;

        return sql.str();
    }

    string row_count_sql(uint matrix_id) const {
        ostringstream sql;
        sql   << "SELECT COUNT(DISTINCT e.i) FROM entries e " << endl
              << "INNER JOIN entries es ON (e.i = es.i) " << endl
              << "WHERE e.matrix_id = " << matrix_id << endl
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ")" << endl
              << " AND  e.type = 'Cell';" << endl;

        return sql.str();
    }

    string row_ids_sql(uint matrix_id) const {
        ostringstream sql;
        sql   << "SELECT DISTINCT e.i FROM entries e " << endl
              << "INNER JOIN entries es ON (e.i = es.i) " << endl
              << "WHERE e.matrix_id = " << matrix_id << endl
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ")" << endl
              << " AND  e.type = 'Cell';" << endl;

        return sql.str();
    }

    string column_ids_sql(uint matrix_id) const {
        ostringstream sql;
        sql   << "SELECT DISTINCT e.j FROM entries e " << endl
              << "INNER JOIN entries es ON (e.i = es.i) " << endl
              << "WHERE e.matrix_id = " << matrix_id << endl
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ")" << endl
              << " AND  e.type = 'Cell';" << endl;

        return sql.str();
    }

    string column_count_sql(uint matrix_id) const {
        ostringstream sql;
        sql   << "SELECT COUNT(DISTINCT e.j) FROM entries e " << endl
              << "INNER JOIN entries es ON (e.i = es.i) " << endl
              << "WHERE e.matrix_id = " << matrix_id << endl
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ")" << endl
              << " AND  e.type = 'Cell';" << endl;

        return sql.str();
    }

    id_set given_ids_;
};

#endif // FUSION_PHENOMATRIX_H_