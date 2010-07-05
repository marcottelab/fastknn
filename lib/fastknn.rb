$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require "distance_matrix.so"

module Fastknn
  VERSION = '0.0.5'

  # Automatically-called function connects to the database. In the future this needs
  # to be revised to take a connection string from Rails.
  def self.connect dbstr = "host=draco.icmb.utexas.edu dbname=crossval_development user=jwoods password=youwish1"
    @@c ||= Fastknn::Connection.new
    @@c.connect(dbstr)
    puts "Connected to database"

    # Create hashes to store matrices, matrix pairs, and distance matrices
    @@source_matrices ||= {}
    @@predict_matrices ||= {}
    @@matrix_pairs ||= {}
    @@distance_matrices ||= {}

    # Keep track of matrices that are cached -- true is cached, non-existent or false is not cached.
    @@cached ||= {}
  end

  # Return a list of cached matrices
  def self.cached
    @@cached.keys.sort
  end

  def self.fetch_source_matrix id, min_genes
    self.mark_as_cached id
    @@source_matrices["#{id}:#{min_genes}"] ||= PhenomatrixBase.new(id, true, min_genes)
  end

  def self.fetch_predict_matrix id, given_id, min_genes
    self.mark_as_cached id
    @@predict_matrices["#{id}:#{given_id}:#{min_genes}"] ||= Phenomatrix.new(id, given_id, min_genes)
  end

  # Cache and return a DistanceMatrix
  def self.fetch_distance_matrix predict_id, source_ids, min_genes
    if source_ids.is_a?(Fixnum)
      source_ids = [source_ids]
    else
      source_ids = source_ids.sort.uniq
    end

    # Cache the matrix pairs
    source_pairs = source_ids.collect { |sid| fetch_matrix_pair(predict_id, sid, min_genes) }

    key = "#{predict_id}:#{source_ids.join(',')}:#{min_genes}"

    self.mark_as_cached id
    @@distance_matrices[key] ||= DistanceMatrix.new(predict_id, source_pairs, min_genes)
  end

  def self.crossvalidate predict_matrix_id, source_matrix_ids, min_genes = 2, distfn = :hypergeometric, classifier_options = {}, dir = "tmp/fastknn"
    opts = {
      :classifier   => :naivebayes,
      :k            => 10,
      :max_distance => 1.0
    }.merge classifier_options

    dm = Fastknn.fetch_distance_matrix(predict_matrix_id, source_matrix_ids, min_genes)
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
  def self.fetch_matrix_pair predict_id, source_id, min_genes
    predict_matrix = fetch_predict_matrix(predict_id, source_id, min_genes)
    source_matrix  = fetch_source_matrix(source_id, min_genes)

    @@matrix_pairs["#{predict_id}:#{source_id}:#{min_genes}"] ||= PhenomatrixPair.new(predict_matrix, source_matrix, min_genes)
  end

  # This allows Rails to lock certain matrices which may be loaded.
  def self.is_cached? matrix_id
    @@cached.has_key?(matrix_id) ? @@cached[matrix_id] : false
  end

protected
  def self.mark_as_cached matrix_id
    @@cached[matrix_id] = true
  end
end

Fastknn.connect
