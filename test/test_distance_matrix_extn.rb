require "rubygems"
require "test_benchmark"

require "test/unit"

require "fastknn"

class TestDistanceMatrixExtn < Test::Unit::TestCase
  def setup
    STDERR.puts "TestDistanceMatrixExtn"
    # Predicting human
    @@d ||= Fastknn.fetch_distance_matrix(185, [3], 2)
    @@d.classifier = {:classifier => :naivebayes, :k => 10, :max_distance => 1}
    @@d.distance_function = :hypergeometric

    d4 ||= Fastknn.fetch_distance_matrix(185, [3], 4)
    d4.classifier = {:classifier => :naivebayes, :k => 10, :max_distance => 1}
    d4.distance_function = :hypergeometric

    
    # Predicting plant
    @@dat ||= Fastknn.fetch_distance_matrix(247, [253,257], 2)
    @@dat.classifier = {:classifier => :naivebayes, :k => 10, :max_distance => 1}
    @@dat.distance_function = :hypergeometric

    @@predict_matrix ||= Fastknn.fetch_predict_matrix(247,247, 2)
    # @@source_matrices ||= @@dat.source_matrix_pairs
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

  def test_min_genes
    d4 = Fastknn.fetch_distance_matrix(185, [3], 4)
    d2 = Fastknn.fetch_distance_matrix(185, [3], 2)
    assert d4.predictable_columns.size == 179
    assert d2.predictable_columns.size == 289

    s4 = Fastknn.fetch_source_matrix(3, 4)
    s2 = Fastknn.fetch_source_matrix(3, 2)
    assert s4.column_ids.size == 3005
    assert s2.column_ids.size == 4056

    j = d2.nearest(3)[0]
    assert j != d4.nearest(3)[0]
    assert s2.has_column?(2950)
    assert !s4.has_column?(2950)
    assert s4.has_column?(4146)
    assert s2.has_column?(4146)
  end

#  def test_predictable_columns
#    s2 = Fastknn.fetch_source_matrix(3, 2)
#    p2 = Fastknn.fetch_matrix_pair(185, 3, 2)
#    d2 = Fastknn.fetch_distance_matrix(185, [3], 2)
#    assert((p2.predictable_columns - s2.column_ids).size == 0)
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
    #assert @@dnearest[1].to_s == "5.97165396725972e-07"
    assert @@dnearest[2] == 3
  end

  def test_predict
    assert @@d.predict(12).size == 15570
    assert @@dat.predict(9831).size == 9388
  end

  def test_push_and_pop_mask
    assert @@dat.predict_matrix_has_column?(9831) == true
    @@dat.push_mask @@first_mask
    assert @@dat.predict_matrix_has_column?(9831) == true
    @@dat.pop_mask
    assert @@dat.predict_matrix_has_column?(9831) == true
    @@dat.push_mask @@second_mask
    assert @@dat.predict_matrix_has_column?(9831) == true
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
