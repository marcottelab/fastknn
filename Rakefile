# File: Rakefile

require 'rubygems'
require 'rdoc/task'
require 'hoe'
require './lib/fastknn'

EXT = "ext/distance_matrix/distance_matrix.so"

Hoe.plugin :git
Hoe.plugin :gemspec
#Hoe.plugin :debugging
#Hoe.plugin :test
#Hoe.plugin :package

h = Hoe.spec 'fastknn' do |p|
  p.developer "John O. Woods", "john.woods@marcottelab.org"
  p.version = '0.0.19'
  p.require_ruby_version '>=1.8'

  p.spec_extras[:extensions] = FileList["ext/**/extconf.rb"]

  p.extra_deps << ["hoe", ">= 2.2.0"]
  p.extra_deps << ["test_benchmark", ">= 0.4.7"]
  p.extra_deps << ["rice", ">= 1.3.2"]
  p.extra_deps << ["rake-compiler", ">= 0.7.0"]

  p.spec_extras[:files] = FileList['Rakefile', 'lib/fastknn.rb', 'test/*.rb', 'ext/**/*.cpp', 'ext/**/*.h'].to_a

  p.clean_globs << EXT << "ext/distance_matrix/*.o" << "ext/distance_matrix/Makefile"

  p.description =<<EOF
Fastknn is a C++-implemented Ruby module for k-nearest neighbors
analyses, providing for multi-stage cross-validation. Includes
a number of classifiers, including hypergeometric, and a bunch
that use TF-IDF such as cosine similarity and Tanimoto co-
efficient.
EOF

  p.need_rdoc = false
end

RDoc::Task.new(:docs) do |rd|
  rd.main = h.readme_file
  rd.options << '-d' if (`which dot` =~ /\/dot/) unless
    ENV['NODOT'] || Hoe::WINDOZE
  rd.rdoc_dir = 'doc'

  rd.rdoc_files.include("lib/**/*.rb")
  rd.rdoc_files += h.spec.extra_rdoc_files
  rd.rdoc_files.reject! {|f| f=="Manifest.txt"}
  title = h.spec.rdoc_options.grep(/^(-t|--title)=?$/).first
  if title then
    rd.options << title

    unless title =~ /\=/ then # for ['-t', 'title here']
    title_index = spec.rdoc_options.index(title)
    rd.options << spec.rdoc_options[title_index + 1]
    end
  else
    title = "#{h.name}-#{h.version} Documentation"
    title = "#{h.rubyforge_name}'s " + title if h.rubyforge_name != h.name
    rd.options << '--title' << title
  end
end

desc 'Publish rdocs with analytics support'
task :publish_docs => [:clean, :docs] do
  ruby %{aggregate_adsense_to_doc.rb}
  path = File.expand_path("~/.rubyforge/user-config.yml")
  config = YAML.load(File.read(path))
  host = "#{config["username"]}@rubyforge.org"

  remote_dir = "/var/www/gforge-projects/#{h.rubyforge_name}/#{h.remote_rdoc_dir
  }"
  local_dir = h.local_rdoc_dir
  Dir.glob(local_dir+"/**/*") {|file|
    sh %{chmod 755 #{file}}
  }
  sh %{rsync #{h.rsync_args} #{local_dir}/ #{host}:#{remote_dir}}
end

require 'rspec/core/rake_task'
namespace :spec do
  desc "Run all specs"
  RSpec::Core::RakeTask.new
  # options in .rspec in package root
end

#namespace :test do
  #Rake::TestTask.new(:rocker) do |t|
  #  t.test_files = FileList['test/test_rocker.rb']
  #  t.warning = true
  #  t.verbose = true
  #end
#end

#Rake::GemPackageTask.new(spec) do |pkg|
#  pkg.need_tar = true
#end

# Rake::ExtensionTask.new('distance_matrix')

file EXT => ["ext/distance_matrix/extconf.rb", "ext/distance_matrix/distance_matrix.cpp"] do
  Dir.chdir "ext/distance_matrix" do
    ruby "extconf.rb"
    sh 'make'
    require "./#{EXT}"
  end
end


#namespace :test do
#  Rake::TestTask.new(:phenomatrix) do |t|
#    t.test_files = FileList['test/test_phenomatrix_extn.rb']
#    t.warning = true
#    t.verbose = true
#  end
#
#  Rake::TestTask.new(:fusion_phenomatrix) do |t|
#    t.test_files = FileList['test/test_fusion_phenomatrix_extn.rb']
#    t.warning = true
#    t.verbose = true
#  end
#
#  Rake::TestTask.new(:phenomatrix_pair) do |t|
#    t.test_files = FileList['test/test_phenomatrix_pair_extn.rb']
#    t.warning = true
#    t.verbose = true
#  end
#
#  Rake::TestTask.new(:distance_matrix) do |t|
#    t.test_files = FileList['test/test_distance_matrix*.rb']
#    t.warning = true
#    t.verbose = true
#  end
#
#  Rake::TestTask.new(:fastknn) do |t|
#    t.test_files = FileList['test/test_fastknn.rb']
#    t.warning = true
#    t.verbose = true
#  end
#end
