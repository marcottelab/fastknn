require 'rubygems'
require 'mkmf-rice'

find_header("phenomatrix.h", "../phenomatrix")
dir_config("distance_matrix")
dir_config("boost")

have_library("stdc++")
have_library("pqxx")
have_library("boost_filesystem")
if RUBY_VERSION =~ /1.9/ then
  $CPPFLAGS += " -DRUBY_19"
end

$CPPFLAGS += " -DRICE"
$objs = ["distance_matrix.o", "hypergeometric.o", "euclidean.o"]

create_makefile('distance_matrix')
