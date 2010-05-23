require "test/unit"

$:.unshift File.dirname(__FILE__) + "/../ext/distance_matrix"
require "distance_matrix.so"

class TestDistanceMatrixExtn < Test::Unit::TestCase
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
    @d ||= DistanceMatrix.new("dbname=crossval_development user=jwoods password=youwish1", 185, [3], "hypergeometric")
    #sz = @d.intersection_size(12, 5143)
    #STDERR.puts "intersection between 12 and 5143 is equal to #{sz}"
    assert @d.intersection_size(12, 5143) == 6
    assert @d.max_intersection_size == 15570
  end

  def test_nearest
    @d ||= DistanceMatrix.new("dbname=crossval_development user=jwoods password=youwish1", 185, [3], "hypergeometric")
    STDERR.puts "Nearest: #{@d.nearest(12)}"
    STDERR.puts "Distance: #{@d.nearest_distance(12)}"
  end

end
