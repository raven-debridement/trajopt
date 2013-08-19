require_relative 'model'
require 'csv'

selectorr = Record.where(version: 19)

#selectorr.each do |record|
#begin
#  record.status = record.result.scan(/^status: (.*)$/).last.last
#  if record.status.strip.size == 0
#    puts record.result.scan(/^status: (.*)$/).inspect
#  end
#rescue
#  puts record.result
#end
#  record.save!
#end

selectorr.where(status: nil).each do |record|
#begin
  record.status = record.result.scan(/^status: (.*)$/).last.last
#rescue
#  puts record.result
#end
  record.save!
end

selectorr.where(run_time: nil).each do |record|
  record.run_time = record.result.scan(/(\d+\.\d+)user/).first.first.to_f
  record.save!
end

selectorr.where(n_merit_increases: nil).each do |record|
  record.n_merit_increases = record.result.scan(/n merit increases: (\d+)/).map{|x| x.last.to_i}.inject(:+)
  record.save!
end

selectorr.where(n_qp_solves: nil).each do |record|
  record.n_qp_solves = record.result.scan(/qp solves: (\d+)/).map{|x| x.last.to_i}.inject(:+)
  record.save!
end

selectorr.where(cost: nil).each do |record|
  record.cost = record.result.scan(/cost values: \((.*), (.*)\)/).last.last.to_f
  record.save!
end

selectorr.where(trust_shrink_ratio: nil).each do |record|
  record.trust_shrink_ratio = record.command.scan(/--trust_shrink_ratio=(\d+\.\d+)/).first.first.to_f
  record.save!
end

selectorr.where(trust_expand_ratio: nil).each do |record|
  record.trust_expand_ratio = record.command.scan(/--trust_expand_ratio=(\d+\.\d+)/).first.first.to_f
  record.save!
end

formulation_map = {
  1 => "stop and turn",
  2 => "constant twist",
}

curvature_constraint_map = {
  1 => "constant curvature",
  2 => "bounded curvature",
}

speed_formulation_map = {
  1 => "constant speed",
  2 => "variable speed",
}

rotation_cost_map = {
  1 => "quadratic rotation cost",
  2 => "L1 rotation cost",
}

selectorr.pluck(:pg_name).uniq.each do |pg_name|

  #[1, 2].product([1, 2], [1, 2], [1, 2]).each do |formulation, curvature_constraint, speed_formulation, rotation_cost|
    selector = selectorr.where(
      pg_name: pg_name,
      #formulation: formulation,
      #curvature_constraint: curvature_constraint,
      #speed_formulation: speed_formulation,
      #rotation_cost: rotation_cost
    )

    converged_count = selector.where(status: "CONVERGED").count
    converged_run_time = selector.where(status: "CONVERGED").pluck(:run_time).inject(:+)
    total_count = selector.count
    total_merit_increases = selector.pluck(:n_merit_increases).inject(:+)
    total_qp_solves = selector.pluck(:n_qp_solves).inject(:+)
    total_cost = selector.pluck(:cost).inject(:+)
    total_run_time = selector.pluck(:run_time).inject(:+)
    #selector.where(:status.ne => "CONVERGED").each do |record|
    #  #puts record.command
    #  #puts record.result.split("\n").last(12)
    #  #puts "STATUS: #{record.status}"
    #end
    #selector.each do |record|
    #  if record.status == "CONVERGED"
    #    converged_count += 1
    #    converged_run_time += record.run_time
    #  #else
    #  #  puts record.command
    #  #  puts record.result.split("\n").last(70)
    #  end
    #  total_run_time += record.run_time
    #  total_merit_increases += record.n_merit_increases
    #  total_qp_solves += record.n_qp_solves
    #  total_cost += record.cost
    #end
    puts "program name: #{pg_name}"
    #puts "formulation: #{formulation_map[formulation]}"
    #puts "curvature_constraint: #{curvature_constraint_map[curvature_constraint]}"
    #puts "speed_formulation: #{speed_formulation_map[speed_formulation]}"
    #puts "rotation_cost: #{rotation_cost_map[rotation_cost]}"
    puts "%converged: #{converged_count.to_f / total_count}"
    puts "average run time for converged: #{converged_run_time / converged_count}"
    puts "average run time overall: #{total_run_time / total_count}"
    puts "average n merit increases: #{total_merit_increases.to_f / total_count}"
    puts "average n qp solves: #{total_qp_solves.to_f / total_count}"
    puts "average cost: #{total_cost.to_f / total_count}"
    puts ""
  #end
end

#CSV.open("results_19.csv", "w") do |csv|
#  csv << %w[program_name formulation curvature_constraint speed_formulation rotation_cost converged_percentage avg_converged_run_time avg_overall_run_time avg_n_merit_increases avg_n_qp_solves avg_cost]
#
#  selectorr.pluck(:pg_name).uniq.each do |pg_name|
#
#    [1, 2].product([1, 2], [1, 2], [1, 2]).each do |formulation, curvature_constraint, speed_formulation, rotation_cost|
#      selector = selectorr.where(
#        pg_name: pg_name,
#        formulation: formulation,
#        curvature_constraint: curvature_constraint,
#        speed_formulation: speed_formulation,
#        rotation_cost: rotation_cost
#      )
#
#      converged_count = selector.where(status: "CONVERGED").count
#      converged_run_time = selector.where(status: "CONVERGED").pluck(:run_time).inject(:+)
#      total_count = selector.count
#      total_merit_increases = selector.pluck(:n_merit_increases).inject(:+)
#      total_qp_solves = selector.pluck(:n_qp_solves).inject(:+)
#      total_cost = selector.pluck(:cost).inject(:+)
#      total_run_time = selector.pluck(:run_time).inject(:+)
#      selector.where(:status.ne => "CONVERGED").each do |record|
#        puts record.command
#        #puts record.result.split("\n").last(12)
#        puts "STATUS: #{record.status}"
#      end
#      #selector.each do |record|
#      #  if record.status == "CONVERGED"
#      #    converged_count += 1
#      #    converged_run_time += record.run_time
#      #  #else
#      #  #  puts record.command
#      #  #  puts record.result.split("\n").last(70)
#      #  end
#      #  total_run_time += record.run_time
#      #  total_merit_increases += record.n_merit_increases
#      #  total_qp_solves += record.n_qp_solves
#      #  total_cost += record.cost
#      #end
#      csv << [pg_name , formulation_map[formulation] , curvature_constraint_map[curvature_constraint] , speed_formulation_map[speed_formulation] , rotation_cost_map[rotation_cost] , converged_count.to_f / total_count , converged_run_time / converged_count , total_run_time / total_count , total_merit_increases.to_f / total_count , total_qp_solves.to_f / total_count , total_cost.to_f / total_count]
#    end
#  end
#end
