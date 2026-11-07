using Rotations
using StaticArrays
using LinearAlgebra
using DifferentialEquations
using DataFrames
using Optim

function phi_p_controller(c::Configuration, shaft_tension::Number, force_h::Number, psi_p::Number)
    # The roll phi tilts the lift out from the centre by T tan(phi). Before the elevation and
    # azimuth rotations the radial unit vector is (0, cos psi_p, sin psi_p), and up across the shaft
    # is (0, 0, -1) (z points down), so tan(phi) = (2 / T)(-F_up sin psi_p + F_h cos psi_p) gives
    # each kite the average force F_up up across the shaft and F_h sideways (the paper's
    # controller, eq. 28-30). Each kite carries the part of the weight across the shaft of itself
    # and of half its tether; the ground carries the other half.
    force_up = (c.m + 0.5 * tether_mass(c.d, c.l)) * c.gravity * cos(c.elev)
    atan(2 / shaft_tension * (-force_up * sin(psi_p) + force_h * cos(psi_p)))
end


function derivative_of_omega_by_psi_p(c::Configuration, wind::Number, omega::Number, psi_p::Number, psi::Number, shaft_tension::Number, mtr::Number, force_h::Number, saver_fun::Function, acc; shaft = shaft_geometry(c, mtr))
    psi_p_list = psi_p .+ range(0, 2 * pi, length = c.n + 1)[1:(end - 1)];
    acc = saver_fun(:psi_p_list, psi_p_list, acc)

    gravity = SVector(0.0, 0.0, c.gravity)
    mass_sum = c.m + tether_mass(c.d, c.l);
    acc = saver_fun(:mass_sum, mass_sum, acc)

    phi_p_list = [phi_p_controller(c, shaft_tension, force_h, p) for p=psi_p_list]
    acc = saver_fun(:phi_p_list, phi_p_list, acc)

    rot_psi_elev = RotZ(psi) * RotY(c.elev)
    rot_psi_elev_psi_p_list = [rot_psi_elev * RotX.(p) for p=psi_p_list]

    c_vec = rot_psi_elev * SVector(1.0, 0.0, 0.0)
    acc = saver_fun(:c_vec, c_vec, acc)
    v_k_vec_list = [rot * SVector(0.0, 0.0, omega * c.radius) for rot=rot_psi_elev_psi_p_list]
    acc = saver_fun(:v_k_vec_list, v_k_vec_list, acc)
    v_a_vec_list = [v .- SVector(wind, 0, 0) for v=v_k_vec_list]
    acc = saver_fun(:v_a_vec_list, v_a_vec_list, acc)
    v_a_list = norm.(v_a_vec_list) .* sign(omega)
    acc = saver_fun(:v_a_list, v_a_list, acc)
    v_a2_list = v_a_list.^2
    r_vec_list = [rot * RotZ(p) * SVector(0.0, 1.0, 0.0) for (rot, p)=zip(rot_psi_elev_psi_p_list, phi_p_list)]
    acc = saver_fun(:r_vec_list, r_vec_list, acc)
    l_n_vec_list = [normalize(cross(r, v)) for (r,v)=zip(r_vec_list, v_a_vec_list)]
    acc = saver_fun(:l_n_vec_list, l_n_vec_list, acc)
    # The lift that holds the shaft tension, less the drag's part along the shaft. The drag follows
    # from C_L, so the lift is found twice: first without the drag, then with the drag of the first.
    lift_for = (drag_along_shaft, l_n) -> (shaft_tension - mass_sum * dot(gravity, c_vec) - drag_along_shaft) / dot(c_vec, l_n)
    c_l_for = (l, v_a2) -> clamp(2 * l / (c.s * c.rho * v_a2), 0.0, c.design_c_l)
    # the kite's polar plus its tether's drag (tether.jl), which depends on the kite's speed and
    # the tension, not on C_L
    c_d_tether_list = [tether_c_d(c, shaft, abs(omega) * c.radius, abs(v_a), shaft_tension) for v_a = v_a_list]
    acc = saver_fun(:c_d_tether_list, c_d_tether_list, acc)
    c_d_of = (c_l, c_d_t) -> calc_c_d(c, c_l, coeffs = c.c_d_fun_coeffs) + c_d_t
    drag_along_shaft = (c_l, c_d_t, v_a, v_a_vec) -> dot(-0.5 * c.rho * v_a * c_d_of(c_l, c_d_t) * c.s * v_a_vec, c_vec)
    c_l_first = [c_l_for(lift_for(0.0, l_n), v_a2) for (l_n, v_a2) = zip(l_n_vec_list, v_a2_list)]
    desired_lift_list = [lift_for(drag_along_shaft(c_l, c_d_t, v_a, v_a_vec), l_n) for (c_l, c_d_t, v_a, v_a_vec, l_n) = zip(c_l_first, c_d_tether_list, v_a_list, v_a_vec_list, l_n_vec_list)]
    acc = saver_fun(:desired_lift_list, desired_lift_list, acc)
    c_l_list = [c_l_for(l, v_a2) for (l, v_a2) = zip(desired_lift_list, v_a2_list)]
    acc = saver_fun(:c_l_list, c_l_list, acc)
    c_d_list = [c_d_of(c_l, c_d_t) for (c_l, c_d_t) = zip(c_l_list, c_d_tether_list)]
    acc = saver_fun(:c_d_list, c_d_list, acc)

    lift_vec_list = [0.5 * c.rho * v_a2 * c_l * c.s * l_n_vec for (v_a2, l_n_vec, c_l) = zip(v_a2_list, l_n_vec_list, c_l_list)]
    acc = saver_fun(:lift_vec_list, lift_vec_list, acc)
    drag_vec_list = [-0.5 * c.rho * v_a * c_d * c.s * v_a_vec for (v_a, v_a_vec, c_d) = zip(v_a_list, v_a_vec_list, c_d_list)]
    acc = saver_fun(:drag_vec_list, drag_vec_list, acc)

    moment_list = [c.radius * dot((l_vec + d_vec + c.m * gravity), normalize(v_k_vec)) for (l_vec, d_vec, v_k_vec) = zip(lift_vec_list, drag_vec_list, v_k_vec_list)]
    acc = saver_fun(:moment_list, moment_list, acc)
    sum_moments = sum(moment_list)
    acc = saver_fun(:sum_moments, sum_moments, acc)
    moment_of_inertia = c.n * c.m * c.radius^2 # tether inertia not accounted for
    acc = saver_fun(:moment_of_inertia, moment_of_inertia, acc)
    # The shaft's moment is Lambda r1 times the tension the kites actually pull, the same tension
    # the power is computed from (it differs from shaft_tension where C_L is clamped).
    shaft_tension_list = [dot(l_vec + d_vec + mass_sum * gravity, c_vec) for (l_vec, d_vec) = zip(lift_vec_list, drag_vec_list)]
    acc = saver_fun(:shaft_tension_list, shaft_tension_list, acc)
    deriv_omega_by_psi_p = 1 / moment_of_inertia / omega * (sum_moments - mtr * c.radius * sum(shaft_tension_list))
    acc = saver_fun(:deriv_omega_by_psi_p, deriv_omega_by_psi_p, acc)
    (deriv_omega_by_psi_p, acc)
end


function solve_sector(c::Configuration, wind::Number, psi::Number, mtr::Number, force_h::Number; step_size::Number = deg2rad(2.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = heuristic_flying_speed(c, wind, psi), shaft_tension::Number = optimal_tension(c, wind, psi, mtr, force_h, step_size = step_size, iterations = iterations, finish_threshold = finish_threshold, speed0 = speed0))
    solver_input = Dict(:config => c
                        , :speed0 => speed0
                        , :wind => wind
                        , :psi => psi
                        , :shaft_tension => shaft_tension
                        , :mtr => mtr
                        , :force_h => force_h
                       )
    psi_p_end = 2 * pi / c.n
    shaft = shaft_geometry(c, mtr)

    deriv = function(du, u, p, psi_p)
        omega = u[1]
        if omega * c.radius < 1.0
            du[1] = 0.0
        else
            du[1] = derivative_of_omega_by_psi_p(c, wind, omega, psi_p, psi, shaft_tension, mtr, force_h, (_, _, _) -> (), (); shaft = shaft)[1]
        end
    end

    solver = function(acc, _)
        if (acc[3]) # done flag
            acc
        else
            omega0 = acc[1]
            u0 = [omega0]
            prob = ODEProblem(deriv, u0, (0, psi_p_end), (solver_input,)) # integrating psi, not time
            # sol = solve(prob, SRIW1(), dt = step_size)
            sol = solve(prob, Euler(), dt = step_size)
            omega1 = sol.u[end][1]
            (omega1, sol, abs(omega0 / omega1 - 1) < finish_threshold)
        end
    end

    (_, result, _) = foldl(solver, 1:max(1, iterations), init = (speed0 / c.radius, (), false))
    result
end


function solve_sector_df(c::Configuration, wind::Number, psi::Number, mtr::Number, force_h::Number; step_size::Number = deg2rad(2.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = heuristic_flying_speed(c, wind, psi), shaft_tension::Number = optimal_tension(c, wind, psi, mtr, force_h, step_size = step_size, iterations = iterations, finish_threshold = finish_threshold, speed0 = speed0), solution::DiffEqBase.AbstractODESolution = solve_sector(c, wind, psi, mtr, force_h, iterations = iterations, step_size = step_size, finish_threshold = finish_threshold, speed0 = speed0, shaft_tension = shaft_tension))
    psi_ps = solution.t
    omegas = [u[1] for u=solution.u]
    shaft = shaft_geometry(c, mtr)
    saved = [derivative_of_omega_by_psi_p(c, wind, omega, psi_p, psi, shaft_tension, mtr, force_h, (k, v, acc) -> (acc[k] = v; acc), Dict{Symbol, Any}(); shaft = shaft)[2] for (psi_p, omega) in zip(psi_ps, omegas)]
    sector_times = diff(psi_ps) ./ omegas[1:(end - 1)]
    t = foldl((acc, _) -> vcat(acc, sector_times .+ maximum(acc)), 1:c.n, init = [0])

    # we want to expand the sector for all 360 degrees
    psi_p_360 = vcat([psi_ps[1:(end - 1)] .+ (i * 2 * pi / c.n) for i=0:(c.n - 1)]...)
    position = [RotZ(psi) * RotY(c.elev) * RotX(psi_p) * SVector(c.l, c.radius, 0.0) for psi_p = psi_p_360]
    saved_value = k -> (tmp = [s[k] for s = saved[1:(end - 1)]]; vcat([tmp for _ = 1:c.n]...))
    saved_value_per_kite = k -> (tmp = [s[k] for s = saved[1:(end - 1)]]; vcat([[row[i] for row=tmp] for i = 1:c.n]...))

    omega = vcat([omegas[1:(end -1)] for _=1:c.n]...)
    mass_sum = saved_value(:mass_sum)
    lift_vec = saved_value_per_kite(:lift_vec_list)
    drag_vec = saved_value_per_kite(:drag_vec_list)
    lift_drag_vec = lift_vec .+ drag_vec
    c_vec = saved_value(:c_vec)
    shaft_tension_actual = [dot(ld + SVector(0.0, 0.0, c.gravity * mass_sum[1]), cv) for (ld, cv) = zip(lift_drag_vec, c_vec)]
    force_h = [dot(ld, RotZ(psi) * SVector(0.0, 1.0, 0.0)) for ld = lift_drag_vec]
    # the aerodynamic force up across the shaft, which carries the weight (z points down)
    up_across_shaft = RotZ(psi) * RotY(c.elev) * SVector(0.0, 0.0, -1.0)
    force_v = [dot(ld, up_across_shaft) for ld = lift_drag_vec]

    # calculate bridle force and angle
    deriv_omega_by_psi_p = saved_value(:deriv_omega_by_psi_p)
    deriv_omega = deriv_omega_by_psi_p .* omega
    v_k_vec = saved_value_per_kite(:v_k_vec_list)
    moment = saved_value_per_kite(:moment_list)
    r_vec = saved_value_per_kite(:r_vec_list)
    tether_force_radial = c.m .* omega .^ 2 .* c.radius .- [dot((lift .+ drag), r_vec) for (lift, drag, r_vec) = zip(lift_vec, drag_vec, r_vec) ] 
    tether_force_travel = c.m .* deriv_omega .* c.radius .- [dot((lift .+ drag), v_k_vec) ./ norm(v_k_vec) - mom / c.radius for (lift, drag, v_k_vec, mom) = zip(lift_vec, drag_vec, v_k_vec, moment) ] 
    bridle_force = [norm([r, t]) for (r, t) = zip(tether_force_radial, tether_force_travel)]
    bridle_angle = [atan(t, r) for (r, t) = zip(tether_force_radial, tether_force_travel)]

    df = DataFrame(:psi_p => psi_p_360
                   , :omega => omega
                   , :power => omega .* mtr .* c.radius .* shaft_tension_actual
                   , :v_k => c.radius .* omega
                   , :t => t[1:(end -1)]
                   , :mass_sum => mass_sum
                   , :phi_p => saved_value_per_kite(:phi_p_list)
                   , :c_vec => c_vec
                   , :v_k_vec => v_k_vec
                   , :v_a_vec => saved_value_per_kite(:v_a_vec_list)
                   , :v_a => saved_value_per_kite(:v_a_list)
                   , :r_vec => r_vec
                   , :l_n_vec => saved_value_per_kite(:l_n_vec_list)
                   , :desired_lift => saved_value_per_kite(:desired_lift_list)
                   , :c_l => saved_value_per_kite(:c_l_list)
                   , :c_d => saved_value_per_kite(:c_d_list)
                   , :c_d_tether => saved_value_per_kite(:c_d_tether_list) # the tether's part of c_d
                   , :lift_vec => lift_vec
                   , :drag_vec => drag_vec
                   , :moment => moment
                   , :sum_moments => saved_value(:sum_moments)
                   , :moment_of_inertia => saved_value(:moment_of_inertia)
                   , :deriv_omega_by_psi_p => deriv_omega_by_psi_p
                   , :deriv_omega => deriv_omega
                   , :position => position
                   , :force_h => force_h
                   , :force_v => force_v
                   , :shaft_tension => shaft_tension_actual
                   , :bridle_force => bridle_force
                   , :bridle_angle => bridle_angle # zero pointing towards center
                  )
    (df, solution.prob.p[1])
end


