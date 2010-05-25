require "rubygems"
require "test_benchmark"

require "test/unit"

$:.unshift File.dirname(__FILE__) + "/../ext/distance_matrix"
require "distance_matrix.so"

class TestDistanceMatrixExtn < Test::Unit::TestCase
  def setup
    @@d ||= Fastknn::DistanceMatrix.new("dbname=crossval_development user=jwoods password=youwish1", 185, [3], "hypergeometric", {:classifier => :naivebayes, :k => 10})
  end

  def test_truth
    assert true
  end

#  def test_distance
#    @d ||= DistanceMatrix.new("dbname=crossval_development user=jwoods password=youwish1", 185, 3, "hypergeometric")
#    dist = @d.distance(7,5,3,10000)
#    # STDERR.puts("Distance is #{dist}")
#    assert 0 < dist
#    assert dist < 1
#  end

  def test_intersection_size
    #sz = @d.intersection_size(12, 5143)
    #STDERR.puts "intersection between 12 and 5143 is equal to #{sz}"
    assert @@d.intersection_size(12, 5143) == 6
  end

  def test_max_intersection_size
    assert @@d.max_intersection_size == 15570
  end

  def test_knearest_size
    @@dnearest ||= @@d.nearest(12)
    knearest = @@d.knearest(12, 10)
    assert knearest.size == 12
  end

  def test_nearest
    @@dnearest ||= @@d.nearest(12)
    assert @@dnearest.first == 1501
    assert @@dnearest[1].to_s == "5.97165396725972e-07"
    assert @@dnearest[2] == 3
  end

  def test_predict
    STDERR.puts "Size of predict(12): #{@@d.predict(12)}"
  end

end
