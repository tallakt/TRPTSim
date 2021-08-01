using Rotations
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using DataFrames

function phi_p_controller(c::Configuration, shaft_tension::Number, force_horizontal::Number, psi_p::Number)
  force_vertical = c.n * (c.m + Tether.mass(c.d, c.l))
  atan(2 / (c.n * shaft_tension) * (force_vertical * sin(psi_p) + force_horizontal * cos(psi_p)))
end


function derivative_of_omega_by_psi_p(c::Configuration, omega::Number, psi_p::Number, psi::Number, shaft_tension::Number, shaft_moment::Number, force_horizontal::Number, saver_fun::Function, acc)
  psi_p_list = psi_p .+ range(0, 2 * pi, length = c.n + 1)[1:(end - 1)];
  acc = saver_fun(:psi_p_list, psi_p_list, acc)

  gravity = SVector(0.0, 0.0, 9.81)
  mass_sum = c.m + Tether.mass(c.d, c.l);
  acc = saver_fun(:mass_sum, mass_sum, acc)

  phi_p_list = [phi_p_controller(c, shaft_tension, force_horizontal, p) for p=psi_p_list]
  acc = saver_fun(:phi_p_list, phi_p_list, acc)

  rot_psi_elev = RotZ(psi) * RotY(c.elev)
  rot_psi_elev_psi_p_list = [rot_psi_elev * RotX.(p) for p=psi_p_list]

  c_vec = rot_psi_elev * SVector(1.0, 0.0, 0.0)
  acc = saver_fun(:c_vec, c_vec, acc)
  v_k_vec_list = [rot * SVector(0.0, 0.0, omega * c.radius) for rot=rot_psi_elev_psi_p_list]
  acc = saver_fun(:v_k_vec_list, v_k_vec_list, acc)
  v_a_vec_list = [v .- SVector(c.w, 0, 0) for v=v_k_vec_list]
  acc = saver_fun(:v_a_vec_list, v_a_vec_list, acc)
  v_a_list = norm.(v_a_vec_list)
  acc = saver_fun(:v_a_list, v_a_list, acc)
  v_a2_list = v_a_list.^2
  r_vec_list = [rot * RotY(p) * SVector(0.0, 1.0, 0.0) for (rot, p)=zip(rot_psi_elev_psi_p_list, phi_p_list)]
  acc = saver_fun(:r_vec_list, r_vec_list, acc)
  l_n_vec_list = [normalize(cross(r, v)) for (r,v)=zip(r_vec_list, v_a_vec_list)]
  acc = saver_fun(:l_n_vec_list, l_n_vec_list, acc)
  lift_list = [(mass_sum * dot(gravity, c_vec) - shaft_tension) / dot(c_vec, l_n) for l_n = l_n_vec_list]
  acc = saver_fun(:lift_list, lift_list, acc)
  c_l_list = [2 * l / (c.s * c.rho * v) for (l, v) = zip(lift_list, v_a2_list)]
  acc = saver_fun(:c_l_list, c_l_list, acc)
  c_d_list = [calc_c_d(c, c_l) for c_l = c_l_list]
  acc = saver_fun(:c_d_list, c_d_list, acc)

  lift_vec_list = [l * l_n for (l, l_n) = zip(lift_list, l_n_vec_list)]
  acc = saver_fun(:lift_vec_list, lift_vec_list, acc)
  drag_vec_list = [-0.5 * c.rho * v_a2 .* c_d .* c.s .* normalize(v_a_vec) for (v_a2, v_a_vec, c_d) = zip(v_a2_list, v_a_vec_list, c_d_list)]
  acc = saver_fun(:drag_vec_list, drag_vec_list, acc)

  moment_list = [c.radius / v_a * dot((l_vec + d_vec + c.m * gravity), v_a_vec) for (v_a, l_vec, d_vec, v_a_vec) = zip(v_a_list, lift_vec_list, drag_vec_list, v_a_vec_list)]
  acc = saver_fun(:moment_list, moment_list, acc)
  sum_moments = sum(moment_list) - shaft_moment
  acc = saver_fun(:sum_moments, sum_moments, acc)
  moment_of_inertia = c.n * c.m * c.radius^2 # tether inertia not accounted for
  acc = saver_fun(:moment_of_inertia, moment_of_inertia, acc)
  deriv_omega_by_psi_p = 1 / moment_of_inertia / omega * sum_moments
  acc = saver_fun(:deriv_omega_by_psi_p, deriv_omega_by_psi_p, acc)
  (deriv_omega_by_psi_p, acc)
end


function solve_sector(c::Configuration, omega0::Number, psi::Number, shaft_tension::Number, shaft_moment::Number, force_horizontal::Number; step_size = deg2rad(1.0))
  solver_input = Dict(:config => c
                      , :omega0 => omega0
                      , :psi => psi
                      , :shaft_tension => shaft_tension
                      , :shaft_moment => shaft_moment
                      , :force_horizontal => force_horizontal
                     )
  psi_p_end = 2 * pi / c.n
  u0 = [omega0]
  deriv = (du, u, p, psi_p) -> du[1] = derivative_of_omega_by_psi_p(c, u[1], psi_p, psi, shaft_tension, shaft_moment, force_horizontal, (_, _, _) -> (), ())[1]
  prob = ODEProblem(deriv, u0, (0, psi_p_end), solver_input) # integrating psi, not time
  solve(prob, Euler(), dt = step_size)
end


function solve_sector_df(c::Configuration, omega0::Number, psi::Number, shaft_tension::Number, shaft_moment::Number, force_horizontal::Number, solution::ODESolution = solve_sector(c, omega0, psi, shaft_tension, shaft_moment, force_horizontal))
  psi_ps = solution.t
  omegas = [u[1] for u=solution.u]
  saved = [derivative_of_omega_by_psi_p(c, omega, psi_p, psi, shaft_tension, shaft_moment, force_horizontal, (k, v, acc) -> (acc[k] = v; acc), Dict{Symbol, Any}())[2] for (psi_p, omega) in zip(psi_ps, omegas)]
  sector_times = diff(psi_ps) ./ omegas[1:(end - 1)]
  t = foldl((acc, _) -> vcat(acc, sector_times .+ maximum(acc)), 1:c.n, init = [0])

  # we want to expand the sector for all 360 degrees
  psi_p_360 = vcat([psi_ps[1:(end - 1)] .+ (i * 2 * pi / c.n) for i=0:(c.n - 1)]...)
  saved_value = k -> (tmp = [s[k] for s = saved[1:(end - 1)]]; vcat([tmp for _ = 1:c.n]...))
  saved_value_per_kite = k -> (tmp = [s[k] for s = saved[1:(end - 1)]]; vcat([[row[i] for row=tmp] for i = 1:c.n]...))
  df = DataFrame(:psi_p => psi_p_360
                 , :omega => vcat([omegas[1:(end -1)] for _=1:c.n]...)
                 , :t => t[1:(end -1)]
                 , :mass_sum => saved_value(:mass_sum)
                 , :phi_p => saved_value_per_kite(:phi_p_list)
                 , :c_vec => saved_value(:c_vec)
                 , :v_k_vec => saved_value_per_kite(:v_k_vec_list)
                 , :v_a_vec => saved_value_per_kite(:v_a_vec_list)
                 , :v_a => saved_value_per_kite(:v_a_list)
                 , :r_vec => saved_value_per_kite(:r_vec_list)
                 , :l_n_vec => saved_value_per_kite(:l_n_vec_list)
                 , :lift => saved_value_per_kite(:lift_list)
                 , :c_l => saved_value_per_kite(:c_l_list)
                 , :c_d => saved_value_per_kite(:c_d_list)
                 , :lift_vec => saved_value_per_kite(:lift_vec_list)
                 , :drag_vec => saved_value_per_kite(:drag_vec_list)
                 , :moment => saved_value_per_kite(:moment_list)
                 , :sum_moments => saved_value(:sum_moments)
                 , :moment_of_inertia => saved_value(:moment_of_inertia)
                 , :deriv_omega_by_psi_p => saved_value(:deriv_omega_by_psi_p)
                )
  (df, solution.prob.p)
end
