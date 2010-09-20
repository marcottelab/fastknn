#!/usr/local/bin/ruby
require "fastknn"

#predict_matrix = Fastknn.fetch_predict_matrix(185,3)
#source_matrix  = Fastknn.fetch_source_matrix(3)
#pair           = Fastknn.fetch_matrix_pair(185,3)
#pair = Fastknn::fetch_matrix_pair(185,3)
d ||= Fastknn.fetch_distance_matrix(1, [1,3,5,7,9,11], 2)
#dm   = Fastknn::fetch_distance_matrix(247,[255],2)
#dm.distance_function = :hyper # clearly wrong
#dm.knearest(9831, 10)

STDERR.puts "Done"
#pair = Fastknn::PhenomatrixPair.new(185,3)
#d ||= Fastknn::DistanceMatrix.new(185, pair)
# p ||= Fastknn.fetch_matrix_pair(185,3)

#d ||= Fastknn.fetch_distance_matrix(185,[3])
#d.classifier = {:classifier => :naivebayes, :k => 10, :max_distance => 1}
#d.distance_function = :hypergeometric
