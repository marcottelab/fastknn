#ifndef DISTANCE_MATRIX_H_
# define DISTANCE_MATRIX_H_

#ifdef RICE
#include <rice/Data_Object.hpp>
#include <rice/Address_Registration_Guard.hpp>
using Rice::Data_Object;
using Rice::Address_Registration_Guard;
#endif

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <map>
#include <set>
#include <algorithm>
#include <string>
using std::string;
using std::set;
using std::endl;
using std::map;
using std::set_intersection;
namespace fs = boost::filesystem;
using fs::ofstream;

typedef std::set<uint> id_set;

#include "fusion_phenomatrix.h"
#include "phenomatrix_pair.h"
#include "cparams.h"
#include "classifier.h"


class Classifier;


class DistanceMatrix {
    friend class Classifier;
    friend class NaiveBayes;
public:

    // This constructor allows a connection to be shared among multiple objects by
    // taking a pointer to an existing one. It is assumed that Ruby won't pass
    // such a pointer, and so we use a set of uints instead of a Rice::Array (which
    // is neccessary for the Ruby interface) as well.
    //
    // In other words, this constructor is exclusively for calling from within
    // a C++ environment of some sort.
    DistanceMatrix(uint, id_set, string, cparams);

    DistanceMatrix(const DistanceMatrix& rhs);

    ~DistanceMatrix();


    // Native C++ cross-validate function. Probably doesn't add much in terms of
    // speed to the native Ruby, but wrote it to simplify debugging.
    void crossvalidate() {
        map<uint, id_set> child_row_sets = predict_matrix_.child_row_ids();

        // Create crossvalidation subdirectories
        prepare_filesystem_for_crossvalidation(child_row_sets.size());

        size_t fold_count = 0;
        for (map<uint, id_set>::const_iterator mt = child_row_sets.begin();
                mt != child_row_sets.end(); ++mt)
        {
            push_mask(mt->second);
            {   // This scope is only for clarity: The predictions here only happen
                // to the masked matrices.
                id_set write_rows(mt->second);

                fs::path dir_name = "predictions" + lexical_cast<string>(fold_count);
                predict_and_write_all_to(dir_name, write_rows);
            }
            pop_mask();

            fold_count++;
        }

    }


    // Removes rows from the matrices on which we're calculating distances.
    void push_mask(id_set mask_rows) {
        for (matrix_list::iterator st = source_matrices.begin(); st != source_matrices.end(); ++st)
            st->push_mask(mask_rows);
    }

    // Restores removed rows from the matrices on which we're calculating distances.
    bool pop_mask() {
        bool res = true;
        for (matrix_list::iterator st = source_matrices.begin(); st != source_matrices.end(); ++st)
            res &= st->pop_mask();
        return res;
    }

    // Find source matrix by columns
    matrix_list::const_iterator find_by_column(uint j) const {
        for (matrix_list::const_iterator dt = source_matrices.begin(); dt != source_matrices.end(); ++dt) {
#ifdef DEBUG_TRACE_INTERSECTION
            cerr << "distance_matrix.h: find_by_column: Checking matrix " << dt->id() << endl;
#endif
            if (dt->source_matrix_has_column(j)) {
#ifdef DEBUG_TRACE_INTERSECTION
                cerr << "\t...found on " << dt->id() << endl;
#endif
                return dt;
            }
        }
        cerr << "distance_matrix.h: Warning: source matrix with column " << j << " was not found." << endl;
        return source_matrices.end(); // not found
    }

    // Determines whether calling observations(j) on the predict matrix is going
    // to throw an exception.
    bool predict_matrix_has_column(uint j) const {
        return predict_matrix_.has_column(j);
    }

    // Find source matrix by matrix ID
    matrix_list::const_iterator find(uint id) const {
        for (matrix_list::const_iterator dt = source_matrices.begin(); dt != source_matrices.end(); ++dt)
            if (dt->id() == id) return dt;
        cerr << "distance_matrix.h: Warning: source matrix with id " << id << " was not found." << endl;
        return source_matrices.end(); // not found
    }


    // Given a row i and a column j, calculate a score (using the classifier function)
    // for the gene's likelihood of being involved in the phenotype.
    pcolumn predict(uint j) const;


    // Predicts for all columns and writes and sorts results
    set<fs::path> predict_and_write_all_to(const fs::path& dir, const id_set& write_rows) const {
        set<fs::path> ret;

        id_set columns = predict_matrix_.column_ids();
        for (id_set::const_iterator jt = columns.begin(); jt != columns.end(); ++jt)
            ret.insert( predict_and_write_to(dir, *jt, write_rows) );

        return ret;
    }

    set<fs::path> predict_and_write_all(id_set write_rows = id_set()) const {
        fs::path dir("predictions");
        return predict_and_write_all_to(dir, write_rows);
    }


    // Write predictions to a file with a specified path. If you don't want to
    // specify a path, use predict_and_write, or set the first argument to "."
    fs::path predict_and_write_to(const fs::path& dir, uint j, const id_set& write_rows) const {
        map<float, id_set> sorted_predictions = predict_rows_and_sort(j, write_rows);

        // Write to a file named by phenotype ID
        fs::path filepath = dir / lexical_cast<string>(j);

        write_predictions(filepath, sorted_predictions);

        return filepath;
    }


    // Predicts for a column, sorts, then returns the filename.
    fs::path predict_and_write(uint j, id_set write_rows = id_set()) const {
        fs::path dir("predictions");
        return predict_and_write_to(dir, j, write_rows);
    }


    // Return the distance, according to our distance function, between the two
    // columns.
    double distance(uint j1, uint j2) const {
        matrix_list::const_iterator j2source = find_by_column(j2);
        if (j2source == source_matrices.end()) {
            string err = "distance_matrix.h: distance(2): column supplied in second argument, " + lexical_cast<string>(j2) + ", does not appear to exist in any source matrix.";
#ifdef RICE
            throw Rice::Exception(rb_eArgError, err.c_str());
#else
            cerr << err << endl;
            throw;
#endif
        }

        return j2source->distance(j1, j2);
    }


    // Note that this function returns only the SINGLE NEAREST -- as in, the first
    // found!
    //
    // If you want multiple nearest, use a different function.
    id_dist_iter nearest(uint j) const {
        // Do a priming read to save an assignment -- particularly since we're
        // likely to have only one source matrix.
        matrix_list::const_iterator pt = source_matrices.begin();
        id_dist_iter min = pt->nearest(j);
        ++pt;

        for (; pt != source_matrices.end(); ++pt) {
            if (!pt->predict_matrix_has_column(j)) continue;

            id_dist_iter min_tmp = pt->nearest(j);
            if (min_tmp.distance < min.distance) {
                min = min_tmp;
                min.matrix_iter = pt; // this keeps us from having to look up each matrix again
            }
        }

        return min;
    }

    // Get the items that are common between j1 in predict matrix and j2 in
    // source matrix
    id_set intersection(uint j1, uint j2) const {
        matrix_list::const_iterator f = find_by_column(j2);
#ifdef DEBUG_TRACE_INTERSECTION
        cerr << "intersection(2): source matrix = " << f->id() << ", j1=" << j1 << ", j2=" << j2 << endl;
#endif
        if (f == source_matrices.end())
            return id_set(); // empty
        
        return f->intersection(j1, j2);
    }

    // Find the k items closest to j
    proximity_queue knearest(uint j, size_t k = 1, double kth_so_far = 100.0) const {
        proximity_queue q;

        for (matrix_list::const_iterator source_matrix_iter = source_matrices.begin(); source_matrix_iter != source_matrices.end(); ++source_matrix_iter) {
            if (!source_matrix_iter->predict_matrix_has_column(j)) continue;
            // Add items on to the queue (q) that are within k (or kth_so_far)
            source_matrix_iter->knearest(q, j, k, kth_so_far, source_matrix_iter);
        }

        // This is the actual return-queue:
        proximity_queue ret;

        // Find the first k items
        while (k > 0 && !q.empty()) {
            kth_so_far = q.top().distance;
            ret.push(q.top());
            q.pop();
            k--;
        }

        // Now include further items equal to kth_so_far
        while (!q.empty() && q.top().distance == kth_so_far) { ret.push(q.top()); q.pop(); }

        return ret;
    }

    // Make a copy and return the prediction matrix as loaded.
    FusionPhenomatrix predict_matrix() const {
        return predict_matrix_;
    }

    id_set predictable_columns() const {
        return predict_matrix_.column_ids();
    }

#ifdef RICE

    // Make a copy and return the list of source matrices (PhenomatrixPairs)
    Rice::Object source_matrix_pairs() const {
        Array ary;
        for (matrix_list::const_iterator i = source_matrices.begin(); i != source_matrices.end(); ++i) {
            ary.push( *i );
        }
        return ary;
    }

    Rice::Object intersection_size(uint i, uint j) const {
        return to_ruby<size_t>(intersection(i, j).size());
    }

    Rice::Object source_matrix_ids() const {
        Rice::Array ids;
        for (matrix_list::const_iterator st = source_matrices.begin(); st != source_matrices.end(); ++st)
            ids.push(to_ruby<uint>(st->id()));
        return ids;
    }

#endif
protected:
    // Create a number of predictions directories in which to store output.
    //
    // Before creation, old directories matching predictions[0-9]+ will be
    // deleted.
    void prepare_filesystem_for_crossvalidation(size_t folds) const {
        // Delete existing directories matching predictions*
        delete_old_crossvalidation_filesystem();

        // Create new directories
        for (size_t f = 0; f < folds; ++f)
            fs::create_directory("predictions" + lexical_cast<string>(f));
    }

    // Delete existing directories matching predictions*
    void delete_old_crossvalidation_filesystem() const {
        fs::path p(".");
        boost::regex e("^predictions[0-9]+$");
        fs::directory_iterator dir_iter(p), dir_end;
        for (; dir_iter != dir_end; ++dir_iter) {
            string filename = dir_iter->filename();
            if (boost::regex_search(filename, e))
                fs::remove_all(filename);
        }
    }

    // Predict for column j and only return the specified set of rows, sorted by score followed by ID
    map<float, id_set> predict_rows_and_sort(uint j, const id_set& rows) const {
        // Simple case -- no rows specified, print all
        if (rows.size() == 0) return predict_and_sort(j);

        // cerr << "predict_rows_and_sort main, rows.size = " << rows.size() << endl;

        pcolumn predictions = predict(j); // first:gene(uint) second:score(float)
        map<float, id_set> sorted_predictions;

        for (id_set::const_iterator i = rows.begin(); i != rows.end(); ++i) {
            pcolumn::const_iterator f = predictions.find(*i);
            // Sort the predictions by adding to a map of sets
            if (f != predictions.end())
                sorted_predictions[f->second].insert(*i);
        }
        return sorted_predictions;
    }

    // Predict for column j and sort the results.
    map<float, id_set> predict_and_sort(uint j) const {
        // cerr << "predict_and_sort main" << endl;

        pcolumn predictions = predict(j); // first:gene(uint) second:score(float)
        map<float, id_set> sorted_predictions;

        // Sort the predictions by adding to a map of sets
        for (pcolumn::const_iterator i = predictions.begin(); i != predictions.end(); ++i)
            sorted_predictions[i->second].insert(i->first);

        return sorted_predictions;
    }
    
    // Output a set of sorted predictions.
    void write_predictions(const fs::path& filename, const map<float, id_set>& sorted_predictions) const {
        ofstream out( filename );

        // Print two header lines to conform with old format
        out << "File generated by fastknn rubygem" << endl;
        out << "gene\tscore" << endl;

        // Print in sorted order: first by score, then by gene if the score is the same
        for (map<float,id_set>::const_reverse_iterator score_it = sorted_predictions.rbegin(); score_it != sorted_predictions.rend(); ++score_it) {
            id_set genes = score_it->second;
            for (id_set::const_iterator gene_it = genes.begin(); gene_it != genes.end(); ++gene_it)
                out << *gene_it << '\t' << score_it->first << endl;
        }

        out.close();
    }


    // Set up the classifier to use for predictions
    void construct_classifier(const cparams&);

    static matrix_list construct_source_matrices(uint predict_matrix_id, const id_set& source_matrix_ids, const string& distfn) {
        matrix_list source_matrices_;
        for (id_set::const_iterator st = source_matrix_ids.begin(); st != source_matrix_ids.end(); ++st) {
            source_matrices_.push_back( PhenomatrixPair(predict_matrix_id, *st, distfn) );
        }
        return source_matrices_;
    }


    // Database-contents-related stuff
    matrix_list source_matrices;
    FusionPhenomatrix predict_matrix_;

    cparams classifier_parameters;
    // Allow different classifiers to be subbed in
    Classifier* classifier;
};


#endif // DISTANCE_MATRIX_H_
