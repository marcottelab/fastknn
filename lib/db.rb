#!/usr/local/bin/ruby
require "fastknn"

#pair = Fastknn::PhenomatrixPair.new(185, 3, "hypergeometric")
#nearest = pair.nearest(12)

#b = Fastknn::Phenomatrix.new(247, 257)
#b.predict(9854)
sids = [257]
pid = 247
#Dir.chdir("tmp/fastknn") do
dm = Fastknn::DistanceMatrix.new pid, sids, "hypergeometric", {:classifier => :naivebayes, :k => 10}
dm.predict(9854)
  #pm = dm.predict_matrix
  #sms = dm.source_matrix_pairs
  #puts dm.predict(12)
  #puts dm.knearest(12,10)
  #puts dm.nearest(12)
#end

