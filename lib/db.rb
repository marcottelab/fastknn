#!/usr/local/bin/ruby
require "fastknn"

#pair = Fastknn::PhenomatrixPair.new(185, 3, "hypergeometric")
#nearest = pair.nearest(12)

sids = [3]
pid = 185
Dir.chdir("tmp/fastknn") do
  dm = Fastknn::DistanceMatrix.new pid, sids, "hypergeometric", {:classifier => :naivebayes, :k => 10}
  #pm = dm.predict_matrix
  #sms = dm.source_matrix_pairs
  #puts dm.predict(12)
  #puts dm.knearest(12,10)
  #puts dm.nearest(12)
end

