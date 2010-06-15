require "rubygems"
require "test_benchmark"

require "test/unit"

# $:.unshift File.dirname(__FILE__) + "/../ext/phenomatrix"
require "fastknn"

class TestPhenomatrixPairExtn < Test::Unit::TestCase
  def setup
    STDERR.puts "TestPhenomatrixPairExtn"
    @@pair ||= Fastknn::PhenomatrixPair.new(185, 3)
    @@pair.distance_function = :hypergeometric
  end

  def test_truth
    assert true
  end

  # NOTE: Distance not guaranteed to be correct. The test is to determine that
  # output is in the correct form and logical.
  def test_push_and_pop_mask
    nearest = @@pair.nearest(12)
    assert nearest[0] == 1501
    assert nearest[1].to_s == "5.97165396725972e-07"
    nearest_before_pop = nearest[1]
    assert nearest[2] == 3

    @@pair.push_mask([675,773,785,2,323,348,351,642,1636,4353,4846])
    nearest = @@pair.nearest(12)
    assert nearest[0] == 2032
    assert nearest[1] > nearest_before_pop
    assert nearest[2] == 3

    @@pair.pop_mask
    nearest = @@pair.nearest(12)
    assert nearest[0] == 1501
    STDERR.puts "nearest to 12 is #{nearest[1].to_s}"
    assert nearest[1] == nearest_before_pop
    assert nearest[2] == 3

    # Shouldn't allow another pop
    assert @@pair.pop_mask == false
  end

end
