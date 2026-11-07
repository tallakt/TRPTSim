using LinearAlgebra
using Roots
using Rotations
using Plots
using LaTeXStrings

# accoring to document chapter about optimal Lambda (moment to tension to radius)
# coordinate system has x-axis horizontal, z-axis along center of shaft and
# y-axis mostly upwards
#
# psi is the heading angle of the kite in the loop, 0 being facing mostly
# upwards, 90 deg being at the top of the loop

function kite_fwd_vec(psi)
  RotZ(psi) * [0.0; 1.0; 0.0]
end


function kite_right_wing_vec(psi)
  RotZ(psi) * [-1.0; 0.0; 0.0]
end


function wind_vec(wind_speed, elevation_angle)
  RotX(-elevation_angle) * [0.0; 0.0; -wind_speed]
end


function apparent_wind(psi, kite_speed, wind_speed, elevation_angle)
  kite_speed * kite_fwd_vec(psi) - wind_vec(wind_speed, elevation_angle)
end

function direction_of_lift_norm_vec(psi, kite_speed, wind_speed, elevation_angle)
  normalize(cross(kite_right_wing_vec(psi), apparent_wind(psi, kite_speed, wind_speed, elevation_angle)))
end


function lift_force_vec(psi, kite_speed, wind_speed, elevation_angle, lift_coefficient, density_of_air, wing_area)
  v_a = norm(apparent_wind(psi, kite_speed, wind_speed, elevation_angle))
  direction = direction_of_lift_norm_vec(psi, kite_speed, wind_speed, elevation_angle)
  0.5 * density_of_air * lift_coefficient * wing_area * v_a^2 * direction
end


function drag_force_vec(psi, kite_speed, wind_speed, elevation_angle, drag_coefficient, density_of_air, wing_area)
  v_a_vec = apparent_wind(psi, kite_speed, wind_speed, elevation_angle)
  v_a = norm(v_a_vec)
  direction = -normalize(v_a_vec)
  0.5 * density_of_air * drag_coefficient * wing_area * v_a^2 * direction
end


function shaft_tension(psi, kite_speed, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
  lift_vec = lift_force_vec(psi, kite_speed, wind_speed, elevation_angle, lift_coefficient, density_of_air, wing_area)
  drag_vec = drag_force_vec(psi, kite_speed, wind_speed, elevation_angle, drag_coefficient, density_of_air, wing_area)
  # never mind harvesting force which is perpendicular to the shaft center line
  direction_of_shaft_vec = [0.0; 0.0; -1.0]
  dot(lift_vec + drag_vec, direction_of_shaft_vec)
end


function lambda(psi, kite_speed, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area, harvesting_force)
  tension = shaft_tension(psi, kite_speed, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
  kite_loop_radius_dummy = 1.0 # just to make the lambda calc easier to understand
  (harvesting_force * kite_loop_radius_dummy) / kite_loop_radius_dummy / tension
end



function find_equilibrium_kite_speed(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area, harvesting_force; initial_speed = wind_speed * lift_coefficient / drag_coefficient * cos(elevation_angle))
  fwd_vec = kite_fwd_vec(psi)
  harvesting_force_vec = -harvesting_force * fwd_vec
  fun = function(kite_speed)
    lift_vec = lift_force_vec(psi, kite_speed, wind_speed, elevation_angle, lift_coefficient, density_of_air, wing_area)
    drag_vec = drag_force_vec(psi, kite_speed, wind_speed, elevation_angle, drag_coefficient, density_of_air, wing_area)
    dot(lift_vec + drag_vec + harvesting_force_vec, fwd_vec)
  end
  try
    find_zero(fun, initial_speed)
  catch
    NaN
  end

end


function maximum_possible_harvesting_force(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
  fun = function(harvesting_force)
    find_equilibrium_kite_speed(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area, harvesting_force)
  end
  max_harvesting_force = find_zero(fun, 0.0)
end


function lambda_power_sweep(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
  max_h_f = maximum_possible_harvesting_force(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
  harvesting_forces = 0:(max_h_f * 0.90)
  kite_speeds = max.(0.0, find_equilibrium_kite_speed.(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area, harvesting_forces))
  shaft_tensions = shaft_tension.(psi, kite_speeds, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
  lambdas = lambda.(psi, kite_speeds, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area, harvesting_forces)
  (lambdas, kite_speeds .* harvesting_forces)
end


println("plot_power_output_vs_lambda_at_different_windspeeds()")
function plot_power_output_vs_lambda_at_different_windspeeds()
  psi = pi / 2
  wind_speeds = 10.0:5.0:30.0
  elevation_angle = pi / 3
  lift_coefficient = 1.5
  drag_coefficient = lift_coefficient / 20.0
  density_of_air = 1.225
  wing_area = 3.0

  p = plot()
  for wind_speed = wind_speeds
    (lambdas, powers) = lambda_power_sweep(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
    plot!(p, lambdas, powers .* 0.001, lab = "$wind_speed m/s", xlabel = L"\Lambda", ylabel = "power [kW]")
  end
  plot!(p, xlims = (0, 0.04), title = L"Power vs $\Lambda$ at wind speeds")
end


println("plot_power_output_vs_lambda_at_different_elevations()")
function plot_power_output_vs_lambda_at_different_elevations()
  psi = pi / 2
  wind_speed = 12.0
  elevation_angles = 0:10:60
  lift_coefficient = 1.5
  drag_coefficient = lift_coefficient / 20.0
  density_of_air = 1.225
  wing_area = 3.0

  p = plot()
  for elevation_angle = elevation_angles
    (lambdas, powers) = lambda_power_sweep(psi, wind_speed, deg2rad(elevation_angle), lift_coefficient, drag_coefficient, density_of_air, wing_area)
    plot!(p, lambdas, powers .* 0.001, lab = "$(elevation_angle)°", xlabel = L"\Lambda", ylabel = "power [kW]")
  end
  plot!(p, xlims = (0, 0.04), title = L"Power vs $\Lambda$ at elevation angles")
end


println("plot_power_output_vs_lambda_at_different_glide_numbers()")
function plot_power_output_vs_lambda_at_different_glide_numbers()
  psi = pi / 2
  wind_speed = 12.0
  elevation_angle = pi / 3
  lift_coefficient = 1.5
  glide_numbers = 5:5:30
  density_of_air = 1.225
  wing_area = 3.0

  p = plot()
  for glide_number = glide_numbers
    (lambdas, powers) = lambda_power_sweep(psi, wind_speed, elevation_angle, lift_coefficient, lift_coefficient / glide_number, density_of_air, wing_area)
    plot!(p, lambdas, powers .* 0.001, lab = "$(glide_number)", xlabel = L"\Lambda", ylabel = "power [kW]")
  end
  plot!(p, xlims = (0, 0.04), title = L"Power vs $\Lambda$ at different $\frac{C_L}{C_D}$")
end



println("plot_power_output_vs_lambda()")
function plot_power_output_vs_lambda()
  psi = pi / 2
  wind_speed = 10.0
  elevation_angle = pi / 3
  lift_coefficient = 1.5
  drag_coefficient = lift_coefficient / 20.0
  density_of_air = 1.225
  wing_area = 3.0

  (lambdas, powers) = lambda_power_sweep(psi, wind_speed, elevation_angle, lift_coefficient, drag_coefficient, density_of_air, wing_area)
  plot(lambdas, powers .* 0.001, lab = "", xlabel = L"\Lambda", ylabel = "power [kW]")
end



nothing
