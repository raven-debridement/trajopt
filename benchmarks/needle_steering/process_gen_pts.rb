cnt = 0
File.open('points_2.txt').read.split("\n").map(&:split).each do |goal_trans_x, goal_trans_y, goal_trans_z|
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


  command = "../../build/bin/check_collision"

  command += options.map {|prefix, val| " #{prefix}#{val}"}.join("")

  result = `#{command}`

  if result.scan(/^status: (.*)$/).first.first == "CONVERGED"
    puts "#{goal_trans_x} #{goal_trans_y} #{goal_trans_z}"
    cnt += 1
  end

  if cnt == 100
    break
  end
end
