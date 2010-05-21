require "test/unit"

$:.unshift File.dirname(__FILE__) + "/../ext/phenomatrix"
require "phenomatrix.so"

class TestPhenomatrixExtn < Test::Unit::TestCase
  def test_truth
    assert true
  end
  
  def test_parent_id
    @p ||= Phenomatrix.new("dbname=crossval_development user=jwoods password=youwish1", 199, true)
    assert @p.parent_id == 193
  end

  def test_parent_id_without_preload
    @p ||= Phenomatrix.new("dbname=crossval_development user=jwoods password=youwish1", 199, false)
    assert @p.parent_id == 193
  end

  def test_root_id
    @p ||= Phenomatrix.new("dbname=crossval_development user=jwoods password=youwish1", 199, true)
    assert @p.root_id == 185
  end

  def test_root_id_without_preload
    @p ||= Phenomatrix.new("dbname=crossval_development user=jwoods password=youwish1", 199, false)
    assert @p.root_id == 185
  end
end
