require 'rubygems'
require 'mkmf-rice'

have_library("stdc++")
have_library("pqxx")
have_library("boost_filesystem")
if RUBY_VERSION =~ /1.9/ then
  $CPPFLAGS += " -DRUBY_19"
end

$CPPFLAGS += " -DRICE"

$objs = ["connection.o"]

create_makefile('connection')
