require 'rubygems'
require 'mkmf-rice'
#require 'mkmf'


have_library("stdc++")
have_library("pqxx")
have_library("boost_filesystem")
# have_library("phenomatrix")
if RUBY_VERSION =~ /1.9/ then
  $CPPFLAGS += " -DRUBY_19"
end

$CPPFLAGS += " -DRICE"
$objs = ["distance_matrix.o", "hypergeometric.o", "euclidean.o", "classifier.o", "phenomatrix_pair.o", "connection.o", "id_dist.o", "naive_bayes.o"]

create_makefile('distance_matrix')
