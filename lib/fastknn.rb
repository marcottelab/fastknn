$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require "lib/phenomatrix.so"

module Fastknn
  VERSION = '0.0.1'

  def self.connect dbstr = "dbname=crossval_development user=jwoods password=youwish1"
    @@c ||= Fastknn::Connection.new
    @@c.connect(dbstr)
    puts "Connected to database"
  end
end

Fastknn.connect