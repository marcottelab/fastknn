require 'rubygems'
require 'mkmf-rice'

dir_config("phenomatrix")
dir_config("boost")

have_library("stdc++")
have_library("pqxx")
have_library("boost")
have_library("boost_filesystem")
if RUBY_VERSION =~ /1.9/ then
  $CPPFLAGS += " -DRUBY_19"
end

$CPPFLAGS += " -DRICE"

create_makefile('phenomatrix')
