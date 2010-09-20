require File.dirname(__FILE__) + '/test_helper.rb'

class TestDistanceMatrixExtn < Test::Unit::TestCase

  def test_min_genes
    source1 = Fastknn.fetch_source_matrix(1,2)
    assert source1.min_genes == 2
    assert !source1.column_ids.include?(284)
    predict1 = Fastknn.fetch_predict_matrix(1, 1, 2)
    assert predict1.min_genes == 2
    assert !predict1.column_ids.include?(284)
    pair1 = Fastknn.fetch_matrix_pair 1, 1, 2
    assert !pair1.predictable_columns.include?(284)
    dm = Fastknn.fetch_distance_matrix 1, [1,3,5,7,9,11], 2
    assert !dm.predictable_columns.include?(284)
  end
end
