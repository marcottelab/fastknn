$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require "distance_matrix.so"

module Fastknn
  VERSION = '0.0.14'

  # Allow Phenomatrix types to be cached as string keys in a hash
  class PhenomatrixBase
    def parent_id?
      self.parent_id.nil? ? nil : self.parent_id
    end

    def inspect additional_h = []
      h = []
      h << ["id", self.id]
      h << ["parent_id", self.parent_id? ? self.parent_id : "nil"]
      h << ["root_id", self.root_id.to_s]
      h << ["row_count", self.row_ids.size]
      h << ["column_count", self.column_ids.size]
      h << ["children_count", self.child_ids.size]
      h << ["min_genes", self.min_genes]
      h.concat additional_h

      "#<#{self.class} " + h.collect { |pair| pair.join(": ") }.join(", ") + ">"
    end

    def to_cache_key
      "#{self.id}:#{self.min_genes}"
    end
  end

  class Phenomatrix
    alias :source_matrix_id :source_id
    alias :given_id :source_id
    def inspect
      super [["source_id", self.source_id]]
    end

    def to_cache_key
      "#{self.id}:#{self.source_matrix_id}:#{self.min_genes}"
    end
  end

  
  class PhenomatrixPair
    alias :predict_matrix_id :predict_id
    alias :source_matrix_id :id
    alias :source_id :id
    def inspect additional_h = []
      h = []
      h << ["id", self.id]
      h << ["predict_matrix_id", self.predict_matrix_id]
      h << ["size", self.size]
      h << ["min_genes", "[#{self.min_genes.join(",")}]"]
      h << ["min_idf", self.min_idf]
      h << ["distance_function", ":#{self.distance_function.to_sym}"]
      h.concat additional_h

      "#<#{self.class} " + h.collect { |pair| pair.join(": ") }.join(", ") + ">"
    end

    def to_cache_key
      "#{self.predict_matrix_id}:#{self.source_matrix_id}:#{self.min_genes}"
    end
  end


  class DistanceMatrix
    alias :predict_matrix_id :id
    def inspect additional_h = []
      dfn_strings = self.distance_functions.collect {|v,d| ":#{d}" }

      h = []
      h << ["id", self.id]
      h << ["source_matrix_ids", "[#{self.source_matrix_ids.join(",")}]"]
      h << ["min_genes", self.min_genes]
      h << ["min_idfs", many_to_one(self.min_idfs.values)]
      h << ["distance_functions", many_to_one(self.distance_functions.values)]
      h.concat additional_h

      "#<#{self.class} " + h.collect { |pair| pair.join(": ") }.join(", ") + ">"
    end

    def to_cache_key
      "#{self.predict_matrix_id}:#{self.source_matrix_ids.sort.join(',')}:#{self.min_genes}"
    end
  protected
    # Take many items that are likely to be the same and describe them as a single
    # item or just 'variable'
    def many_to_one arr
      if arr.uniq.size == 1
        arr.first
      else
        '(variable)'
      end
    end
  end

  # Automatically-called function connects to the database. In the future this needs
  # to be revised to take a connection string from Rails.
  def self.connect dbstr = "host=localhost dbname=crossval_development user=jwoods password=youwish1"
    @@c ||= Fastknn::Connection.new
    @@c.connect(dbstr)
    puts "Connected to database"

    # Create hashes to store matrices, matrix pairs, and distance matrices
    @@source_matrices ||= {}
    @@predict_matrices ||= {}
    @@matrix_pairs ||= {}
    @@distance_matrices ||= {}

    # Keep track of matrices that are cached -- true is cached, non-existent or false is not cached.
    @@cached ||= Hash.new{ |h,k| h[k] = [] }
  end

  # Return a list of cached matrices
  def self.cached
    @@cached.keys.sort
  end

  def self.fetch_source_matrix id, min_genes
    key = "#{id}:#{min_genes}"

    @@source_matrices[key] ||= PhenomatrixBase.new(id, true, min_genes)
    item = @@source_matrices[key]

    self.mark_as_cached id, item

    item
  end

  def self.fetch_predict_matrix id, given_id, min_genes
    key = "#{id}:#{given_id}:#{min_genes}"

    @@predict_matrices[key] ||= Phenomatrix.new(id, given_id, min_genes)
    item = @@predict_matrices[key]

    self.mark_as_cached id, item
    self.mark_as_cached given_id, item

    item
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

    @@distance_matrices[key] ||= DistanceMatrix.new(predict_id, source_pairs, min_genes)
    item = @@distance_matrices[key]

    self.mark_as_cached predict_id, item
    source_ids.each { |given_id|  self.mark_as_cached(given_id, item) }

    item
  end

  def self.crossvalidate predict_matrix_id, source_matrix_ids, min_genes = 2, distfn = :hypergeometric, min_idf = 0.0, classifier_options = {}, dir = "tmp/fastknn"
    opts = {
      :classifier   => :naivebayes,
      :k            => 10,
      :max_distance => 1.0
    }.merge classifier_options

    dm = Fastknn.fetch_distance_matrix(predict_matrix_id, source_matrix_ids, min_genes)
    dm.classifier        = opts
    dm.distance_function = distfn
    dm.min_idf           = min_idf

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

    key                   = "#{predict_id}:#{source_id}:#{min_genes}"
    @@matrix_pairs[key] ||= PhenomatrixPair.new(predict_matrix, source_matrix, min_genes)
    
    item                  = @@matrix_pairs[key]
    Fastknn.mark_as_cached predict_id, item
    
    item
  end

  # This allows Rails to lock certain matrices which may be loaded.
  def self.is_cached? matrix_id
    @@cached.has_key?(matrix_id) ? @@cached[matrix_id] : false
  end

  def self.uncache matrix_id
    @@cached[matrix_id].each do |cache_item|
      case cache_item.class.to_s
      when "Fastknn::PhenomatrixBase" then @@source_matrices.delete(cache_item.to_cache_key)
      when "Fastknn::Phenomatrix"     then @@predict_matrices.delete(cache_item.to_cache_key)
      when "Fastknn::PhenomatrixPair" then @@matrix_pairs.delete(cache_item.to_cache_key)
      when "Fastknn::DistanceMatrix"  then @@distance_matrices.delete(cache_item.to_cache_key)
      else
        STDERR.puts "Error uncaching, unknown class #{cache_item.class}"
        return false
      end
    end

    @@cached.delete(matrix_id)
    true
  end

protected
  def self.mark_as_cached id, cache_item
    @@cached[id] << cache_item
  end
end

Fastknn.connect
