require 'mongoid'

Mongoid.load!("./mongoid.yml", :development)

class Record
  include Mongoid::Document

  field :status
  field :run_time
  field :n_merit_increases
  field :n_qp_solves
  field :cost
  field :trust_shrink_ratio
  field :trust_expand_ratio
end
