require "rubygems"
require "test_benchmark"

require "test/unit"

#$:.unshift File.dirname(__FILE__) + "/../ext/phenomatrix"
require "fastknn"

class TestPhenomatrixExtn < Test::Unit::TestCase
  def setup
    @@p ||= Fastknn::Phenomatrix.new(185, 3)
    @@pb ||= Fastknn::PhenomatrixBase.new(185)
    @@pp ||= Fastknn::Phenomatrix.new(185, 185)
  end

  def test_truth
    assert true
  end
  
  def test_row_count
    assert @@pp.row_count == 16648
    assert @@pb.row_count == 16648
    assert @@p.row_count == 15570
  end

  def test_parent_and_root_id
    @@p193 ||= Fastknn::PhenomatrixBase.new(199)
    assert @@p193.parent_id == 193
    assert @@p.root_id == 185
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

  def test_child_ids
    child_ids = @@p.child_ids
    assert child_ids.size == 5
    assert child_ids.include?(193)
    assert child_ids.include?(194)
    assert child_ids.include?(195)
    assert child_ids.include?(196)
    assert child_ids.include?(197)
  end

  def test_child_row_ids
    child_row_ids = @@pp.child_row_ids
    assert child_row_ids.size == 5
    assert child_row_ids[193].size == 3330
    assert child_row_ids[194].size == 3330
    assert child_row_ids[195].size == 3330
    assert child_row_ids[196].size == 3329
    assert child_row_ids[197].size == 3329
  end

end
