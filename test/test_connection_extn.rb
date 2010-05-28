require "test/unit"

$:.unshift File.dirname(__FILE__) + "/../ext/connection"
require "connection.so"

class TestConnectionExtn < Test::Unit::TestCase
  def test_truth
    assert true
  end
end