# File: Rakefile

require 'rubygems'
require 'rake/extensiontask'
require 'rake/testtask'

spec = Gem::Specification.new do |s|
  s.name = "Fastknn"
  s.platform = Gem::Platform::RUBY
  s.extensions = FileList["ext/**/extconf.rb"]
  s.version = "0.0.1"
end

Rake::GemPackageTask.new(spec) do |pkg|
end

Rake::ExtensionTask.new('phenomatrix')
Rake::ExtensionTask.new('distance_matrix')

namespace :test do
Rake::TestTask.new(:all) do |t|
  t.test_files = FileList['test/test*.rb']
  t.warning = true
  t.verbose = true
end

Rake::TestTask.new(:phenomatrix) do |t|
  t.test_files = FileList['test/test_phenomatrix*.rb']
  t.warning = true
  t.verbose = true
end

Rake::TestTask.new(:distance_matrix) do |t|
  t.test_files = FileList['test/test_distance_matrix*.rb']
  t.warning = true
  t.verbose = true
end
end
# Rake::ExtensionTask.new('classifier')
