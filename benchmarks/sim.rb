require 'parallel'
require_relative 'model'

sigmas = (0.25..3).step(0.25).to_a
noise_levels = [1]

options = [
  { open_loop: 1, use_lqr: 1 },
  { open_loop: 1, use_lqr: 0 },
  { open_loop: 0, use_lqr: 0 },
]

100.times do
  seed = rand(100000000)
  while Record.where(seed: seed).count > 0
    seed = rand(100000000)
  end
  Parallel.each(sigmas.product(noise_levels, options), :in_threads => 8) do |sigma, noise_level, options|
    open_loop = options[:open_loop]
    use_lqr = options[:use_lqr]
    command = "time ../build/bin/four_links_robot_sim --initial_controls_path=../build/controls/control_#{sigma} --sigma_pts_scale=#{sigma} --noise_level=#{noise_level} --seed=#{seed} --open_loop=#{open_loop} --use_lqr=#{use_lqr} 2>&1"
    result = `#{command}`
    Record.create! sigma: sigma, noise_level: noise_level, seed: seed, result: result, open_loop: !open_loop.zero?, use_lqr: !use_lqr.zero?, command: command
  end
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
