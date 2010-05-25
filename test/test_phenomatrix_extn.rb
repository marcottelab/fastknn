require "rubygems"
require "test_benchmark"

require "test/unit"

$:.unshift File.dirname(__FILE__) + "/../ext/phenomatrix"
require "phenomatrix.so"

class TestPhenomatrixExtn < Test::Unit::TestCase
  def setup
    @@p ||= Fastknn::Phenomatrix.new("dbname=crossval_development user=jwoods password=youwish1", 199)
  end

  def test_truth
    assert true
  end
  
  def test_parent_id
    assert @@p.parent_id == 193
  end

  def test_root_id
    assert @@p.root_id == 185
  end

end
