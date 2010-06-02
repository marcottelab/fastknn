require "rubygems"
require "test_benchmark"

require "test/unit"

require "fastknn"

class TestDistanceMatrixExtn < Test::Unit::TestCase
  def setup
    @@d ||= Fastknn::DistanceMatrix.new(185, [3], "hypergeometric", {:classifier => :naivebayes, :k => 10})
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
    assert @@d.predict(12).size == 16648
  end

  def test_push_and_pop_mask
    @@d.push_mask([675,773,785,2,323,348,351,642,1636,4353,4846])
    @@d.pop_mask

    @@d.push_mask([1,9, 12, 13, 15, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32, 34, 35, 36])
    @@d.pop_mask
  end

#  def test_predict_and_write
#    Dir.chdir("tmp/fastknn") do
#      @@d.predict_and_write(1)
#      assert `cat 1 |wc -l` =~ /^16650$/
#      `rm 1`
#
#      @@d.predict_and_write_all([1,9,12,13,15,19,20,22,23])
#      assert `cat 4 |wc -l` =~ /^11$/
#      `rm *`
#
#      @@d.predict_and_write_all
#      assert `cat 5 |wc -l` =~ /^16650$/
#      `rm *`
#    end
#  end

  def test_crossvalidation
    Fastknn.crossvalidate 185, 3, "hypergeometric", {:classifier => :naivebayes, :k => 10}
  end

end
