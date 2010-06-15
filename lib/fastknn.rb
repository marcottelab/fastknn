$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require "distance_matrix.so"

module Fastknn
  VERSION = '0.0.4'

  # Automatically-called function connects to the database. In the future this needs
  # to be revised to take a connection string from Rails.
  def self.connect dbstr = "dbname=crossval_development user=jwoods password=youwish1"
    @@c ||= Fastknn::Connection.new
    @@c.connect(dbstr)
    puts "Connected to database"

    # Create hashes to store matrices, matrix pairs, and distance matrices
    @@source_matrices ||= {}
    @@predict_matrices ||= {}
    @@matrix_pairs ||= {}
    @@distance_matrices ||= {}
  end

  def self.fetch_source_matrix id
    @@source_matrices[id] ||= PhenomatrixBase.new(id)
  end

  def self.fetch_predict_matrix id, given_id
    @@predict_matrices["#{id}:#{given_id}"] ||= Phenomatrix.new(id, given_id)
  end

  # Cache and return a DistanceMatrix
  def self.fetch_distance_matrix predict_id, source_ids
    if source_ids.is_a?(Fixnum)
      source_ids = [source_ids]
    else
      source_ids = source_ids.sort.uniq
    end

    # Cache the matrix pairs
    source_pairs = source_ids.collect { |sid| fetch_matrix_pair(predict_id, sid) }

    key = "#{predict_id}:#{source_ids.join(',')}"

    @@distance_matrices[key] ||= DistanceMatrix.new(predict_id, source_pairs)
  end

  def self.crossvalidate predict_matrix_id, source_matrix_ids, distfn = :hypergeometric, classifier_options = {}, dir = "tmp/fastknn"
    opts = {
      :classifier   => :naivebayes,
      :k            => 10,
      :max_distance => 1.0
    }.merge classifier_options

    dm = Fastknn.fetch_distance_matrix(predict_matrix_id, source_matrix_ids)
    dm.classifier        = opts
    dm.distance_function = distfn

    puts "Current dir = #{Dir.pwd}"

    # Make sure directory exists
    unless File.exist?(dir) && File.directory?(dir)
      Dir.new dir
    end

    # Change to that directory and run the cross-validation function
    Dir.chdir(dir) do   
      dm.crossvalidate
    end
  end

#protected
  # Cache and return a PhenomatrixPair. This is protected because we don't want
  # the user pushing or popping masks.
  def self.fetch_matrix_pair predict_id, source_id
    predict_matrix = fetch_predict_matrix(predict_id, source_id)
    source_matrix  = fetch_source_matrix(source_id)

    @@matrix_pairs["#{predict_id}:#{source_id}"] ||= PhenomatrixPair.new(predict_matrix, source_matrix)
  end

end

Fastknn.connect
