cnt = 0
File.open('points_100000.txt').read.split("\n").map(&:split).each do |goal_trans_x, goal_trans_y, goal_trans_z|
  goal_rot_x = 0
  goal_rot_y = 0
  goal_rot_z = 0

  options = []
  options << ['--g ', goal_trans_x]
  options << ['--g ', goal_trans_y]
  options << ['--g ', goal_trans_z]
  options << ['--g ', goal_rot_x]
  options << ['--g ', goal_rot_y]
  options << ['--g ', goal_rot_z]

  options << ['--T=', 1]
  options << ['--continuous_collision=', 0]
  options << ['--control_constraints=', 0]
  options << ['--collision_dist_pen=', 0.05]
  options << ['--r_min=', 0.1]


  command = "../../build/bin/needle_steering"

  command += options.map {|prefix, val| " #{prefix}#{val}"}.join("")

  result = `#{command}`

  if result.scan(/^status: (.*)$/).first.first == "CONVERGED"
    puts "#{goal_trans_x} #{goal_trans_y} #{goal_trans_z}"
    cnt += 1
  end

  if cnt == 10000
    break
  end
end
