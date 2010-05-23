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

Rake::TestTask.new(:test_all) do |t|
  t.test_files = FileList['test/test*.rb']
  t.warning = true
  t.verbose = false
end
# Rake::ExtensionTask.new('classifier')
