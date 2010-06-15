require "rubygems"
require "test_benchmark"
require "set"

require "test/unit"

# $:.unshift File.dirname(__FILE__) + "/../ext/phenomatrix"
require "fastknn"

class TestFusionPhenomatrixExtn < Test::Unit::TestCase
  def setup
    STDERR.puts "TestPhenomatrixExtn"
    @@fp ||= Fastknn::FusionPhenomatrix.new(247, [253,257])
    @@p253 ||= Fastknn::Phenomatrix.new(247, 253, 2)
    @@p257 ||= Fastknn::Phenomatrix.new(247, 257, 2)
    # @@pair ||= Fastknn::PhenomatrixPair.new(247, [253,257], "hypergeometric")
  end

  def test_truth
    assert true
  end

  def test_row_count_triangularity
    pair_row_set = Set.new(@@p253.row_ids).merge(@@p257.row_ids)
    fusion_row_set = Set.new(@@fp.row_ids)
    STDERR.puts("p253: #{@@p253.row_ids.size}")
    STDERR.puts("p257: #{@@p257.row_ids.size}")
    STDERR.puts("union: #{pair_row_set.size}")
    STDERR.puts("fusion: #{fusion_row_set.size}")

    assert fusion_row_set == pair_row_set
  end

end
