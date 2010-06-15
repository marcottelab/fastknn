#!/usr/local/bin/ruby
require "fastknn"

#predict_matrix = Fastknn.fetch_predict_matrix(185,3)
#source_matrix  = Fastknn.fetch_source_matrix(3)
#pair           = Fastknn.fetch_matrix_pair(185,3)
pair = Fastknn::fetch_matrix_pair(185,3)
dm   = Fastknn::fetch_distance_matrix(185,3)
STDERR.puts "Done"
#pair = Fastknn::PhenomatrixPair.new(185,3)
#d ||= Fastknn::DistanceMatrix.new(185, pair)
# p ||= Fastknn.fetch_matrix_pair(185,3)

#d ||= Fastknn.fetch_distance_matrix(185,[3])
#d.classifier = {:classifier => :naivebayes, :k => 10, :max_distance => 1}
#d.distance_function = :hypergeometric
