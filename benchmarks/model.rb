require 'mongoid'

Mongoid.load!("./mongoid.yml", :development)

class Record
  include Mongoid::Document
end
