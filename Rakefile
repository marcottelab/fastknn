# File: Rakefile

require 'rubygems'
require 'rake/extensiontask'
require 'rake/testtask'
require 'hoe'

PKG = "fastknn"
PKG_VERSION = [0,0,4]
AUTHOR = "John O. Woods, Marcotte Lab"
EMAIL = "john.woods@marcottelab.org"
HOMEPAGE = "http://github.com/MarcotteLabGit/fastknn"

Hoe.plugin :git
Hoe.plugin :debugging
Hoe.plugin :test

Hoe.spec PKG do |p|
  p.developer AUTHOR, EMAIL

  p.readme_file    = "README.rdoc"

  p.extra_deps << ["hoe", ">= 2.2.0"]
  p.extra_deps << ["test_benchmark", ">= 0.4.7"]
  p.extra_deps << ["rice", ">= 1.3.2"]
end

#spec = Gem::Specification.new do |s|
#  s.platform = Gem::Platform::RUBY
#  s.extensions = FileList["ext/**/extconf.rb"]
#  s.summary = "Faster k-nearest neighbors analysis for crossval"
#  s.name = PKG
#  s.author = AUTHOR
#  s.email = EMAIL
#  s.homepage = HOMEPAGE
#  s.version = PKG_VERSION.join('.')
#  s.requirements << 'libpqxx3'
#  s.requirements << 'rice-1.3.2'
#  s.require_path = 'lib'
#  # s.autorequire = 'rake'
#  s.files = FileList['Rakefile', 'lib/fastknn.rb', 'test/*.rb', 'ext/**/*.cpp', 'ext/**/*.h'].to_a
#  s.description = <<EOF
#Fastknn is a C++-implemented Ruby module for k-nearest neighbors
#analyses, providing for multi-stage cross-validation.
#EOF
#end
#spec.add_development_dependency('test_benchmark')
#
#Rake::GemPackageTask.new(spec) do |pkg|
#  pkg.need_tar = true
#end
#
Rake::ExtensionTask.new('distance_matrix')


namespace :test do
  Rake::TestTask.new(:phenomatrix) do |t|
    t.test_files = FileList['test/test_phenomatrix_extn.rb']
    t.warning = true
    t.verbose = true
  end

  Rake::TestTask.new(:fusion_phenomatrix) do |t|
    t.test_files = FileList['test/test_fusion_phenomatrix_extn.rb']
    t.warning = true
    t.verbose = true
  end

  Rake::TestTask.new(:phenomatrix_pair) do |t|
    t.test_files = FileList['test/test_phenomatrix_pair_extn.rb']
    t.warning = true
    t.verbose = true
  end

  Rake::TestTask.new(:distance_matrix) do |t|
    t.test_files = FileList['test/test_distance_matrix*.rb']
    t.warning = true
    t.verbose = true
  end
end
