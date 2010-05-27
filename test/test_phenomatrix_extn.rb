require "rubygems"
require "test_benchmark"

require "test/unit"

$:.unshift File.dirname(__FILE__) + "/../ext/phenomatrix"
require "phenomatrix.so"

class TestPhenomatrixExtn < Test::Unit::TestCase
  def setup
    @@p ||= Fastknn::Phenomatrix.new("dbname=crossval_development user=jwoods password=youwish1", 185, 3)
    @@pb ||= Fastknn::PhenomatrixBase.new("dbname=crossval_development user=jwoods password=youwish1", 185)
    @@pp ||= Fastknn::Phenomatrix.new("dbname=crossval_development user=jwoods password=youwish1", 185, 185)
  end

  def test_truth
    assert true
  end
  
  def test_parent_and_root_id
    @@p193 ||= Fastknn::PhenomatrixBase.new("dbname=crossval_development user=jwoods password=youwish1", 199)
    assert @@p193.parent_id == 193
    assert @@p.root_id == 185
  end

  def test_successful_database_load
    assert @@p.row_count > 0
    assert @@pp.row_count > 0
    assert @@pb.row_count > 0
  end

  # If we use dest = source, row count should be exactly the same. So should everything.
  def test_phenomatrix_base_and_phenomatrix_equivalence
    assert @@pb.row_count == 16648
    assert @@pb.row_count == @@pp.row_count
  end

  # If we do not use dest = source, we expect the two-species phenomatrix to have
  # fewer rows than the one-species phenomatrix
  def test_phenomatrix_base_and_phenomatrix_inequivalence
    assert @@p.row_count < @@pb.row_count
    assert @@p.row_count == 15570
  end

end
