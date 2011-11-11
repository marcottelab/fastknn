# -*- encoding: utf-8 -*-

Gem::Specification.new do |s|
  s.name = "fastknn"
  s.version = "0.0.19.20111111163231"

  s.required_rubygems_version = Gem::Requirement.new(">= 0") if s.respond_to? :required_rubygems_version=
  s.authors = ["John O. Woods"]
  s.date = "2011-11-11"
  s.description = "Fastknn is a C++-implemented Ruby module for k-nearest neighbors\nanalyses, providing for multi-stage cross-validation. Includes\na number of classifiers, including hypergeometric, and a bunch\nthat use TF-IDF such as cosine similarity and Tanimoto co-\nefficient.\n"
  s.email = ["john.woods@marcottelab.org"]
  s.extensions = ["ext/distance_matrix/extconf.rb"]
  s.extra_rdoc_files = ["History.txt", "Manifest.txt", "PostInstall.txt"]
  s.files = ["Rakefile", "lib/fastknn.rb", "test/test_phenomatrix_extn.rb", "test/test_fastknn.rb", "test/test_distance_matrix_extn.rb", "test/test_phenomatrix_pair_extn.rb", "test/test_helper.rb", "test/test_fusion_phenomatrix_extn.rb", "test/test_connection_extn.rb", "ext/distance_matrix/rice_connection.cpp", "ext/distance_matrix/simple_classifier.cpp", "ext/distance_matrix/phenomatrix_base.cpp", "ext/distance_matrix/similarity.cpp", "ext/distance_matrix/ruby_conversions.cpp", "ext/distance_matrix/hypergeometric.cpp", "ext/distance_matrix/fusion_phenomatrix.cpp", "ext/distance_matrix/id_dist.cpp", "ext/distance_matrix/connection.cpp", "ext/distance_matrix/rice_distance_matrix.cpp", "ext/distance_matrix/phenomatrix_pair.cpp", "ext/distance_matrix/phenomatrix.cpp", "ext/distance_matrix/classifier.cpp", "ext/distance_matrix/distance_matrix.cpp", "ext/distance_matrix/average_classifier.cpp", "ext/distance_matrix/rice_phenomatrix.cpp", "ext/distance_matrix/naive_bayes.cpp", "ext/distance_matrix/similarity.h", "ext/distance_matrix/typedefs.h", "ext/distance_matrix/id_dist.h", "ext/distance_matrix/classifier.h", "ext/distance_matrix/connection.h", "ext/distance_matrix/phenomatrix_pair.h", "ext/distance_matrix/phenomatrix.h", "ext/distance_matrix/params.h", "ext/distance_matrix/hypergeometric.h", "ext/distance_matrix/fusion_phenomatrix.h", "ext/distance_matrix/distance_matrix.h", "History.txt", "Manifest.txt", "PostInstall.txt", "ext/distance_matrix/extconf.rb", ".gemtest"]
  s.homepage = "http://github.com/marcottelab/fastknn"
  s.rdoc_options = ["--main", "README.txt"]
  s.require_paths = ["lib"]
  s.required_ruby_version = Gem::Requirement.new(">= 1.8")
  s.rubyforge_project = "fastknn"
  s.rubygems_version = "1.8.10"
  s.summary = "FIX (describe your package)"
  s.test_files = ["test/test_phenomatrix_extn.rb", "test/test_fastknn.rb", "test/test_distance_matrix_extn.rb", "test/test_phenomatrix_pair_extn.rb", "test/test_helper.rb", "test/test_fusion_phenomatrix_extn.rb", "test/test_connection_extn.rb"]

  if s.respond_to? :specification_version then
    s.specification_version = 3

    if Gem::Version.new(Gem::VERSION) >= Gem::Version.new('1.2.0') then
      s.add_runtime_dependency(%q<hoe>, [">= 2.2.0"])
      s.add_runtime_dependency(%q<test_benchmark>, [">= 0.4.7"])
      s.add_runtime_dependency(%q<rice>, [">= 1.3.2"])
      s.add_runtime_dependency(%q<rake-compiler>, [">= 0.7.0"])
      s.add_development_dependency(%q<hoe>, ["~> 2.12"])
    else
      s.add_dependency(%q<hoe>, [">= 2.2.0"])
      s.add_dependency(%q<test_benchmark>, [">= 0.4.7"])
      s.add_dependency(%q<rice>, [">= 1.3.2"])
      s.add_dependency(%q<rake-compiler>, [">= 0.7.0"])
      s.add_dependency(%q<hoe>, ["~> 2.12"])
    end
  else
    s.add_dependency(%q<hoe>, [">= 2.2.0"])
    s.add_dependency(%q<test_benchmark>, [">= 0.4.7"])
    s.add_dependency(%q<rice>, [">= 1.3.2"])
    s.add_dependency(%q<rake-compiler>, [">= 0.7.0"])
    s.add_dependency(%q<hoe>, ["~> 2.12"])
  end
end
