require 'parallel'
require_relative 'model'

def rand_range(min, max)
  rand * (max - min) + min
end

T = 15
cnt = 0

r_min = 4.5

File.open('new_points_10000.txt').read.split("\n").map(&:split).each do |goal_trans_x, goal_trans_y, goal_trans_z|
  cnt += 1
  if cnt > 250
    break
  end
  puts "simulation No.#{cnt}"

  %w[needle_steering_grid needle_steering_single_stage].each do |pg_name|

    #dis = rand_range(5, 10)
    #theta = rand_range(Math::PI/8, Math::PI/4)
    
    goal_rot_x = 0
    goal_rot_y = 0
    goal_rot_z = 0

    #r_min = rand_range(0.4, 0.5)# * dis
    options = []
    options << ['--r_min=', r_min]

    #if pg_name == "needle_steering"
      options << ['--improve_ratio_threshold=',  0.1]
      options << ['--trust_shrink_ratio=', 0.9]
      options << ['--trust_expand_ratio=', 1.3]
    #end
    
    options << ['--T=', T]
    #options << ['--s ', start_trans_x]
    #options << ['--s ', start_trans_y]
    #options << ['--s ', start_trans_z]
    #options << ['--s ', start_rot_x]
    #options << ['--s ', start_rot_y]
    #options << ['--s ', start_rot_z]
    options << ['--g ', goal_trans_x]
    options << ['--g ', goal_trans_y]
    options << ['--g ', goal_trans_z]
    options << ['--g ', goal_rot_x]
    options << ['--g ', goal_rot_y]
    options << ['--g ', goal_rot_z]

    curvature_formulation = 1
    
    [1, 2].product([1, 2], [1, 2], [1, 2]).each do |formulation, curvature_constraint, speed_formulation, rotation_cost|
    #[1].product([1], [1], [1], [0.1, 1, 10, 100, 1000]).each do |formulation, curvature_constraint, speed_formulation, rotation_cost, coeff_orientation_error|
    #[1].product([1], [1], [1], [1]).each do |formulation, curvature_constraint, speed_formulation, rotation_cost, coeff_orientation_error|
      exp_options = options.clone

      exp_options << ['--formulation=', formulation]
      exp_options << ['--curvature_constraint=', curvature_constraint]
      exp_options << ['--curvature_formulation=', curvature_formulation]
      exp_options << ['--speed_formulation=', speed_formulation]
      if speed_formulation == 2
        exp_options << ['--use_speed_deviation_cost=', 1]
      end
      exp_options << ['--rotation_cost=', rotation_cost]
      exp_options << ['--collision_dist_pen=', 0.05]
      exp_options << ['--method=', 1]
      exp_options << ['--plot_final_result=', 0]
      exp_options << ['--plotting=', 0]

      command = "time ../../build/bin/#{pg_name}"
      command += exp_options.map {|prefix, val| " #{prefix}#{val}"}.join("")
      command += " 2>&1"

      puts command

      #exit(1)
      result = `#{command}`

      #puts result

      Record.create! formulation: formulation,
                     curvature_constraint: curvature_constraint,
                     curvature_formulation: curvature_formulation,
                     speed_formulation: speed_formulation,
                     rotation_cost: rotation_cost,
                     method: 1,
                     result: result,
                     command: command,
                     pg_name: pg_name,
                     version: 29,
                     description: "single stage vs grid with different configurations"

    end
  end
end
