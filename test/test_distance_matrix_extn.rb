require "rubygems"
require "test_benchmark"

require "test/unit"

require "fastknn"

class TestDistanceMatrixExtn < Test::Unit::TestCase
  def setup
    # Predicting human
    @@d ||= Fastknn::DistanceMatrix.new(185, [3], "hypergeometric", {:classifier => :naivebayes, :k => 10})
    
    # Predicting plant
    @@dat ||= Fastknn::DistanceMatrix.new(247, [253,257], "hypergeometric", {:classifier => :naivebayes, :k => 10})

    @@predict_matrix ||= Fastknn::Phenomatrix.new(247,247)
    @@source_matrices ||= @@dat.source_matrix_pairs
    @@masks ||= @@predict_matrix.child_row_ids
    @@first_mask ||= @@masks[264]
    @@second_mask ||= @@masks[259]
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
    assert @@dat.predict_matrix_has_column?(9831) == true
    assert @@dat.predict_matrix_has_column?(9854) == false
    @@dat.push_mask @@first_mask
    assert @@dat.predict_matrix_has_column?(9831) == true
    assert @@dat.predict_matrix_has_column?(9854) == false
    @@dat.pop_mask
    assert @@dat.predict_matrix_has_column?(9831) == true
    assert @@dat.predict_matrix_has_column?(9854) == false
    @@dat.push_mask @@second_mask
    assert @@dat.predict_matrix_has_column?(9831) == true
    assert @@dat.predict_matrix_has_column?(9854) == false
    @@dat.pop_mask
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

#  def test_crossvalidation
#    Fastknn.crossvalidate 185, 3, "hypergeometric", {:classifier => :naivebayes, :k => 10}
#  end

end
