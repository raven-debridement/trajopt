require_relative 'model.rb'

puts "#{Record.count} record(s) in total"

Record.where(:in_collision => nil).each do |record|
  if record.result =~ /NOT IN COLLISION/
    record.update_attributes(:in_collision => false)
  else
    record.update_attributes(:in_collision => true)
  end
end

Record.where(:distance_to_goal => nil).each do |record|
begin
  record.update_attributes(:distance_to_goal => /distance to goal: (.*)$/.match(record.result)[1].to_f)
rescue
  puts record.result
end
end

puts "collision statistics (w/ replanning):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.25..3).step(0.25)
              .map{|sigma| [sigma, Record.where(sigma: sigma, open_loop: false, use_lqr: false)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}

puts "collision statistics (open loop execution):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.25..3).step(0.25)
              .map{|sigma| [sigma, Record.where(sigma: sigma, open_loop: true, use_lqr: false)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}

puts "collision statistics (lqr):"
puts "sigma\tpercentage_in_collision\taverage_distance_to_goal"
puts (0.25..3).step(0.25)
              .map{|sigma| [sigma, Record.where(sigma: sigma, open_loop: true, use_lqr: true)]}
              .map{|sigma, selector| [sigma, selector.count, selector.where(in_collision:true).count, selector.pluck(:distance_to_goal).inject(&:+)]}
              .map{|sigma, count, collision_count, sum_distance| "#{sigma}\t#{collision_count.to_f/count}\t#{sum_distance.to_f/count}"}
