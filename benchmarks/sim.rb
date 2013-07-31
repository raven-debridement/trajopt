require 'parallel'
require_relative 'model'

sigmas = [0]#(0.5..4).step(0.5).to_a
noise_levels = [1]

options = [
  { open_loop: 1, use_lqr: 1, initial_solve: 0 },
  #{ open_loop: 1, use_lqr: 0, initial_solve: 0 },
  { open_loop: 0, use_lqr: 0, initial_solve: 0 },
]

#seeds = []

seeds = [99854284, 23217141, 19434567, 42815875, 4009578, 12810802, 46972549, 92877962, 16105370, 43507157, 98487417, 92980143, 17660850, 65900921, 25376427, 14495540, 91034473, 26149965, 52868298, 82131720, 77472759, 93377342, 15044138, 38934700, 24089759, 94564002, 76766353, 11880620, 23174001, 50215280, 278460, 78434736, 34141769, 43985220, 74749021, 41017832, 22233645, 99058424, 80059478, 65672992, 5270312, 12301743, 37686166, 48194559, 42388776, 1709759, 44066982, 24462569, 59263333, 58928555, 53986640, 89322204, 31712734, 60716006, 55891718, 60397700, 33292074, 24625600, 2852263, 5808473, 19001639, 10443045, 48125694, 26404100, 14008686, 23688635, 68243107, 90379126, 27654828, 68319716, 89790489, 98629283, 79043711, 18978928, 29214086, 29627575, 34305293, 8106975, 60988870, 35331951, 55247156, 22515301, 21697692, 64826541, 9209814, 40108585, 93189270, 91912758, 43741763, 9600903, 80452796, 62597200, 4401855, 93915294, 9542028, 36412447, 88995204, 91704711, 39055649, 66153025]
#100.times do
#  seed = rand(100000000)
#  while Record.where(seed: seed).count > 0 or seeds.include? seed
#    seed = rand(100000000)
#  end
#  seeds << seed
#end
#
#puts seeds.inspect
  

#cnt = 0
#Record.each do |record|
#  if not record.result =~ /distance to goal:/
#    puts "running"
#    cnt += 1
#    record.update_attributes result: `#{record.command}`
#  end
#end
#
#puts cnt

Parallel.each(seeds.product(sigmas, noise_levels, options), :in_threads => 8) do |seed, sigma, noise_level, options|
  open_loop = options[:open_loop]
  use_lqr = options[:use_lqr]
  initial_solve = options[:initial_solve]
  command = "time ../build/bin/four_links_robot_sim --initial_controls_path=../build/controls/state_space --sigma_pts_scale=#{sigma} --noise_level=#{noise_level} --seed=#{seed} --open_loop=#{open_loop} --use_lqr=#{use_lqr} --initial_solve=#{initial_solve} 2>&1"
  result = `#{command}`
  #puts result
  Record.create! sigma: sigma, noise_level: noise_level, seed: seed, result: result, type: :state_space, open_loop: !open_loop.zero?, use_lqr: !use_lqr.zero?, initial_solve: !initial_solve.zero?, command: command
end

#Record.pluck(:seed).uniq.each do |seed|
#  Parallel.each(Record.where(:open_loop => false, :use_lqr => false, :seed => seed), :in_threads => 8) do |record|
#    sigma = record.sigma
#    noise_level = record.noise_level
#    seed = record.seed
#    result = `time ../build/bin/four_links_robot --sigma_pts_scale=#{sigma} --noise_level=#{noise_level} --seed=#{seed} --open_loop=1 --use_lqr=1 2>&1`
#    Record.create! sigma: sigma, open_loop: true, use_lqr: true, noise_level: noise_level, seed: seed, result: result
#  end
#end
