#ifndef FUSION_PHENOMATRIX_H_
# define FUSION_PHENOMATRIX_H_

#include "phenomatrix.h"

#include <sstream>
#include <string>
#include <algorithm>

using std::string;
using std::ostringstream;
using std::insert_iterator;

#ifndef MATRIX_LIST_DEFINED
# define MATRIX_LIST_DEFINED
# include "phenomatrix_pair.h"
typedef std::list<PhenomatrixPair>          matrix_list;
#endif

id_set extract_matrix_ids(const matrix_list& matrices);
id_set extract_row_ids(const matrix_list& matrix_pairs);
id_set extract_column_ids(const matrix_list& matrix_pairs);
id_set union_ids(id_set a, id_set b);
id_set intersect_ids(id_set a, id_set b);


class FusionPhenomatrix : public PhenomatrixBase {
public:
    FusionPhenomatrix(uint id, id_set given_ids, size_t min_genes = 2)
    : PhenomatrixBase(id, true, min_genes),
            given_ids_(given_ids)
    {
        if (given_ids.size() > 1 || *(given_ids.begin()) != id)
            inherit_construct(min_genes);
    }

    FusionPhenomatrix(uint id, const matrix_list& given_matrices);

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
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ");" << endl;

        return sql.str();
    }

    string row_ids_sql(uint matrix_id) const {
        ostringstream sql;
        sql   << "SELECT DISTINCT e.i FROM entries e " << endl
              << "INNER JOIN entries es ON (e.i = es.i) " << endl
              << "WHERE e.matrix_id = " << matrix_id << endl
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ") ORDER BY e.i;" << endl;

        return sql.str();
    }

    string column_ids_sql(uint matrix_id) const {
        ostringstream sql;
        sql   << "SELECT DISTINCT e.j FROM entries e " << endl
              << "INNER JOIN entries es ON (e.i = es.i) " << endl
              << "WHERE e.matrix_id = " << matrix_id << endl
              << " AND es.matrix_id IN (" << join(given_ids_, ",") << ")" << endl
              << " AND  e.type = 'Cell' order by e.j;" << endl;

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

