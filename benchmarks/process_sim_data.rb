require_relative 'model.rb'

puts "#{Record.count} record(s) in total"

Record.where(:in_collision => nil).each do |record|
#Record.each do |record|
  if record.result =~ /NOT IN COLLISION/
    record.update_attributes(:in_collision => false)
  elsif record.result =~ /IN COLLISION/
    record.update_attributes(:in_collision => true)
  else
    puts record.result
  end
end

Record.where(:distance_to_goal => nil).each do |record|
#Record.each do |record|
begin
  record.update_attributes(:distance_to_goal => /distance to goal: (.*)$/.match(record.result)[1].to_f)
rescue
  puts record.result
end
end

selector = Record.where(open_loop: false, use_lqr: false, type: :belief_space)
puts "collision statistics (belief space w/ replanning, #{selector.count} in total):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.5..4).step(0.5)
              .map{|sigma| [sigma, selector.where(sigma: sigma)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}

selector = Record.where(open_loop: true, use_lqr: false, type: :belief_space)
puts "collision statistics (belief space open loop execution, #{selector.count} in total)):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.5..4).step(0.5)
              .map{|sigma| [sigma, selector.where(sigma: sigma)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}

selector = Record.where(open_loop: true, use_lqr: true, type: :belief_space)
puts "collision statistics (belief space lqr, #{selector.count} in total)):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.5..4).step(0.5)
              .map{|sigma| [sigma, selector.where(sigma: sigma)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}

selector = Record.where(open_loop: false, use_lqr: false, type: :state_space)
puts "collision statistics (state space w/ replanning, #{selector.count} in total):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.0..0.0).step(0.5)
              .map{|sigma| [sigma, selector.where(sigma: sigma)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}

selector = Record.where(open_loop: true, use_lqr: true, type: :state_space)
puts "collision statistics (state space lqr, #{selector.count} in total)):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.0..0.0).step(0.5)
              .map{|sigma| [sigma, selector.where(sigma: sigma)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}
