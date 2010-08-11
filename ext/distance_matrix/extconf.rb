require 'rubygems'
require 'mkmf-rice'
#require 'mkmf'


have_library("stdc++")
have_library("pqxx")
have_library("boost_filesystem")
have_library("boost_regex")
# have_library("phenomatrix")
if RUBY_VERSION =~ /1.9/ then
  $CPPFLAGS += " -DRUBY_19"
end

$CPPFLAGS += " -DRICE -DNDEBUG"
$CXXFLAGS += " -O3 -DNDEBUG"
#$CPPFLAGS += " -DRICE -g3 -O0"
#$CXXFLAGS += " -O0 -g3"

$objs = [
  "distance_matrix.o",
  "hypergeometric.o",
  "classifier.o",
  "simple_classifier.o",
  "average_classifier.o",
  "naive_bayes.o",
  "phenomatrix.o",
  "phenomatrix_pair.o",
  "connection.o",
  "id_dist.o",
  "fusion_phenomatrix.o",
  "similarity.o"]

create_makefile('distance_matrix')
