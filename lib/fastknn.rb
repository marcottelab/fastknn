$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require "lib/distance_matrix.so"

module Fastknn
  VERSION = '0.0.1'

  # Automatically-called function connects to the database. In the future this needs
  # to be revised to take a connection string from Rails.
  def self.connect dbstr = "dbname=crossval_development user=jwoods password=youwish1"
    @@c ||= Fastknn::Connection.new
    @@c.connect(dbstr)
    puts "Connected to database"
  end

  def self.crossvalidate predict_matrix_id, source_matrix_ids, distfn = "hypergeometric", classifier_options = {}
    opts = {
      :classifier => :naivebayes,
      :k          => 10
    }.merge classifier_options

    source_matrix_ids = [source_matrix_ids] if source_matrix_ids.is_a?(Fixnum)

    predict_matrix = predict_matrix_id
    predict_matrix = Matrix.find(predict_matrix_id) unless predict_matrix.is_a?(Matrix)

    dm = DistanceMatrix.new(predict_matrix_id, source_matrix_ids, distfn, opts)

    counter = 1

    list_of_row_sets_for(predict_matrix).each do |row_set|
      # Tell it which rows to ignore in the distance calculations.
      dm.push_mask row_set

      STDERR.puts("Predicting #{counter}")

      # Tell it which rows to predict and write to files.
      dm.predict_and_write_all row_set

      # Now remove the mask we added before.
      dm.pop_mask

      counter += 1
    end
  end


  # Given a matrix, find the children and get the rows associated. 
  def self.list_of_row_sets_for predict_matrix
    # Get list of row sets
    row_sets = []
    predict_matrix.children.each do |child|
      row_sets << child.rows
    end

    row_sets
  end

end

Fastknn.connect