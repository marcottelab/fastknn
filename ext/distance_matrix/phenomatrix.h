#ifndef PHENOMATRIX_H_
# define PHENOMATRIX_H_

#include <string>
#include <boost/lexical_cast.hpp>
using std::string;
using boost::lexical_cast;

#include "phenomatrix_base.h"

// This derivative class of PhenomatrixBase will typically represent the prediction
// species in a pair of species. It can be used for either! Just set the two ids
// to the same value.
//
// As described in PhenomatrixBase's comment section, a prediction matrix will
// have more rows than a source matrix. This Phenomatrix takes that into account
// when loading the prediction matrix and only loads the part which is relevant
// to the source matrix.
class Phenomatrix : public PhenomatrixBase {
public:
    // id = predict species, given_id = source species
    Phenomatrix(uint id, uint given_id, size_t min_genes = 2)
    : PhenomatrixBase(id, id == given_id, min_genes),
      given_id_(given_id)
    {
        // this is the part of the construction that is different between the
        // base and derived class.
        if (id != given_id) inherit_construct(min_genes);
    }

    Phenomatrix(const Phenomatrix& rhs)
    : PhenomatrixBase(rhs), given_id_(rhs.given_id_)
    { }

    // Copy and mask
    Phenomatrix(const Phenomatrix& rhs, const id_set& mask_rows)
    : PhenomatrixBase(rhs, mask_rows), given_id_(rhs.given_id_)
    { }

    ~Phenomatrix() { }

    // Get the source matrix id
    uint source_id() const { return given_id_; }
protected:
    // OVERRIDE SQL FROM PHENOMATRIX_BASE.
    string load_matrix_sql(uint matrix_id) const {
        return string("SELECT DISTINCT e1.i,e1.j FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + " \n" +
               " AND  e1.type = 'Cell' ORDER BY e1.j, e1.i;";
    }

    string row_count_sql(uint matrix_id) const {
        return string("SELECT COUNT(DISTINCT e1.i) FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + ";";
    }

    string row_ids_sql(uint matrix_id) const {
        return string("SELECT DISTINCT e1.i FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + " ORDER BY e1.i;";
    }

    string column_ids_sql(uint matrix_id) const {
        return string("SELECT DISTINCT e1.j FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + " ORDER BY e1.j;";
    }

    string column_count_sql(uint matrix_id) const {
        return string("SELECT COUNT(DISTINCT e1.j) FROM entries e1 \n") +
               "INNER JOIN entries e2 ON (e1.i = e2.i) \n" +
               "WHERE e1.matrix_id = " + lexical_cast<string>(matrix_id) + " \n" +
               " AND  e2.matrix_id = " + lexical_cast<string>(given_id_) + ";";
    }

    uint given_id_;
};

#endif