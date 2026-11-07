using TRPTSim
using Test
using Statistics
using LinearAlgebra

@testset "Very simplified kite model" begin
  glide = 10.0
  c = config(n = 3
             , elev = 0.0
             , l = 30.0
             , c_d_tether = 0.0
             , s = 1.0
             , d = 0.01                   # thick but no drag
             , design_c_l = 2.0
             , c_d_fun_coeffs = [0.001, 1 / glide]
             , m = 2.0                     # light, but not so light that the Euler steps get unstable
             , radius = 15.0

           )
  w = 12.0
  psi = 0.0
  m_t_factor = 0.0
  force_h = 0.0
  tension_factor = 0.5
  # want to fly with half possible tension
  tension = tension_factor * (0.5 * c.rho * c.design_c_l * c.s * glide ^ 2 * w ^ 2)
  (df, solve_input) = solve_sector_df(c, w, psi, m_t_factor, force_h, shaft_tension = tension)

  @test mean(df.v_a) ≈ norm([1, glide]) * w             atol = 2.0
  @test mean(df.shaft_tension) ≈ tension                atol = 100.0
  @test mean(signal_sum_of_kites(df.shaft_tension, c.n)) ≈ tension * c.n    atol = 300.0
  @test mean(df.c_l) ≈ c.design_c_l * tension_factor    atol = 0.02
end


@testset "Very simplified kite model and production" begin
  glide = 10.0
  c = config(n = 3
             , elev = 0.0
             , l = 30.0
             , c_d_tether = 0.0
             , s = 1.0
             , d = 0.01                   # thick but no drag
             , design_c_l = 2.0
             , c_d_fun_coeffs = [0.001, 1 / glide]
             , m = 2.0                     # light, but not so light that the Euler steps get unstable
             , radius = 15.0

           )
  w = 12.0
  psi = 0.0
  m_t_factor = 0.02 # the kite then flies near w / (Lambda + 1 / glide), where it can hold the tension
  force_h = 0.0
  tension_factor = 0.5
  # want to fly with half possible tension
  tension = tension_factor * (0.5 * c.rho * c.design_c_l * c.s * glide ^ 2 * w ^ 2)
  (df, solve_input) = solve_sector_df(c, w, psi, m_t_factor, force_h, shaft_tension = tension)

  @test maximum(df.c_l) < c.design_c_l
  # m_t_factor is Lambda, the moment per tension per radius
  @test mean(df.moment) ≈ tension .* m_t_factor .* c.radius                     rtol = 0.01
  @test mean(df.power) ≈ tension .* m_t_factor .* c.radius .* mean(df.omega)    rtol = 0.01
end

@testset "Tether and tether drag tests" begin
  glide = 10.0
  c = config(n = 3
             , elev = 0.0
             , l = 30.0
             , c_d_tether = 1.0
             , s = 1.0
             , d = 0.01                   # thick but no drag
             , design_c_l = 2.0
             , c_d_fun_coeffs = [0.001, 1 / glide]
             , m = 2.0                     # light, but not so light that the Euler steps get unstable
             , radius = 15.0

           )
  w = 12.0
  psi = 0.0
  m_t_factor = 0.001
  force_h = 0.0
  tension_factor = 0.5
  # want to fly with half possible tension
  tension = tension_factor * (0.5 * c.rho * c.design_c_l * c.s * glide ^ 2 * w ^ 2)
  (df, solve_input) = solve_sector_df(c, w, psi, m_t_factor, force_h, shaft_tension = tension)

  # the tether adds C_t d l / (4 S) to the drag coefficient, whatever the kite's C_L
  c_d_tether = c.c_d_tether * c.d * c.l / 4 / c.s
  for c_l in (0.5, 1.0, 2.0)
    @test calc_c_d(c, c_l) - calc_c_d(c, c_l, coeffs = c.c_d_fun_coeffs) ≈ c_d_tether
  end
  # flying free at elevation 0, the kite flies at its glide number, v_a = w sqrt(1 + G^2), with the
  # drag coefficient it flew with (the tether's part depends on the speed and the tension)
  @test df.c_d ≈ [calc_c_d(c, cl, coeffs = c.c_d_fun_coeffs) for cl = df.c_l] .+ df.c_d_tether
  glide_flown = mean(df.c_l) / mean(df.c_d)
  @test mean(df.v_a) ≈ w * sqrt(1 + glide_flown^2)    rtol = 0.03
end



@testset "Test elevation and phi" begin
  w = 12.0
  m_t_factor = 0.025
  force_h = 0.0
  c = config(eijkelhof(), gravity = 0.0)
  elev_psi = [(0.0, 0.0), (30.0, 0.0), (60.0, 0.0), (0.0, 30.0), (0.0, 60.0), (30.0, 30.0)]
  powers = [solve_sector_df(config(c, elev = deg2rad(elv)), w, deg2rad(psi), m_t_factor, force_h) |> get_avg_power for (elv, psi) = elev_psi]
  # without gravity only the angle between the wind and the shaft counts
  angle_30_30 = acos(cos(deg2rad(30.0)) * cos(deg2rad(30.0)))
  power_30_30 = solve_sector_df(config(c, elev = angle_30_30), w, 0.0, m_t_factor, force_h) |> get_avg_power

  @test powers[1] > powers[2] > powers[3]
  @test powers[4] ≈ powers[2]   rtol = 0.01
  @test powers[5] ≈ powers[3]   rtol = 0.01
  @test powers[6] ≈ power_30_30 rtol = 0.01
end


@testset "Test scaling and effect of mass" begin
  w = 9.0
  m_t_factor = 0.025
  force_h = 0.0
  psi = 0.0
  c1 = eijkelhof()
  c2 = config(scale_to_area(c1, c1.s / 2), m = c1.m / 2)
  c3 = scale_to_area(c1, c1.s / 2)
  power1 = solve_sector_df(c1, w, psi, m_t_factor, force_h) |> get_avg_power
  power2 = solve_sector_df(c2, w, psi, m_t_factor, force_h) |> get_avg_power
  power3 = solve_sector_df(c3, w, psi, m_t_factor, force_h) |> get_avg_power

  @test c3.m < c2.m # mass should scale more than x^2
  @test power1 ≈ 2 * power2     rtol = 0.05 # power goes as the area
  @test power3 ≈ power2         rtol = 0.100
  @test power3 > power2 # holding up the weight across the shaft costs power
end


@testset "Test horizontal and vertical force generation" begin
  w = 9.0
  m_t_factor = 0.025
  psi = 0.0
  tension = 300_000.0
  c = eijkelhof()
  (sol1, _) = solve_sector_df(c, w, psi, m_t_factor, 0.0, shaft_tension = tension)
  (sol2, _) = solve_sector_df(c, w, psi, m_t_factor, 1000.0, shaft_tension = tension)
  (sol3, _) = solve_sector_df(c, w, psi, m_t_factor, 2000.0, shaft_tension = tension)
  (sol4, _) = solve_sector_df(c, w, psi, m_t_factor, -2000.0, shaft_tension = tension)

  # each kite holds up its weight across the shaft and half its tether's; the kite's drag, which
  # the controller leaves out, takes a little of it
  force_up = (c.m + tether_mass(c.d, c.l) / 2) * c.gravity * cos(c.elev)
  @test mean(sol1.force_v) ≈ force_up rtol = 0.25
  @test mean(signal_sum_of_kites(sol1.force_v, c.n)) ≈ c.n * mean(sol1.force_v) rtol = 0.001
  @test mean(sol1.force_h) ≈ 0.0      atol = 600.0
  @test mean(sol2.force_h) ≈ 1000.0   atol = 700.0
  @test mean(sol3.force_h) ≈ 2000.0   atol = 700.0
  @test mean(sol4.force_h) ≈ -mean(sol3.force_h) atol = 1500.0
end


@testset "Roll controller holds the shaft up" begin
  # a kite with almost no drag, so that the lift alone makes the force across the shaft
  c = config(eijkelhof(), c_d_tether = 0.0, c_d_fun_coeffs = [0.0002, 0.002])
  (df, _) = solve_sector_df(c, 9.0, 0.0, 0.025, 0.0, shaft_tension = 300_000.0)
  force_up = (c.m + tether_mass(c.d, c.l) / 2) * c.gravity * cos(c.elev)
  @test mean(df.force_v) > 0 # up, not down
  # the controller sets the roll from the paper's averages, open loop, which hold to about 10 %
  @test mean(df.force_v) ≈ force_up rtol = 0.15
  @test maximum(df.phi_p) > 0 && minimum(df.phi_p) < 0
  # at the bottom of the loop (psi_p = 90 degrees, z points down) the kite rolls its lift inwards, up
  i_bottom = argmin(abs.(df.psi_p .- pi / 2))
  @test df.force_v[i_bottom] > 0
end


@testset "The shaft is braked by the tension it carries" begin
  c = eijkelhof()
  tension = 300_000.0
  # a tight periodic orbit, so that the ring's kinetic energy is the same at the start and the end
  (df, input) = solve_sector_df(c, 9.0, 0.0, 0.025, 0.0, shaft_tension = tension, step_size = deg2rad(0.5),
                                finish_threshold = 1e-9, iterations = 5000)
  @test maximum(df.c_l) < c.design_c_l
  # with the drag along the shaft in the lift, the kites pull the tension asked for
  @test mean(df.shaft_tension) ≈ tension rtol = 0.001
  # over a loop the generator gets the work the air does on the kites (gravity and the speed
  # change of the ring give nothing back over a period); averages over time, dt = dpsi / omega
  time_mean = x -> sum(x ./ df.omega) / sum(1 ./ df.omega)
  aero_work = [dot(l + d, v) for (l, d, v) = zip(df.lift_vec, df.drag_vec, df.v_k_vec)]
  @test time_mean(df.power) ≈ time_mean(aero_work) rtol = 0.001
end




@testset "A tension is found in strong wind" begin
  # a kite with a high glide number, whose heuristic tension is far above what its tether carries
  c = config(october_kite(polar = :profile), c_d_tether = 0.0)
  t = optimal_tension(c, 20.0, 0.0, 0.025, 0.0)
  @test !isnan(t)
  @test t <= 1.01 * tether_strength(c.d) / c.safety_factor
end


@testset "The October kite's polar" begin
  # the physical polar (the default) is the October polar plus the induced drag of the wing,
  # span 26 m on 54 m^2 (aspect ratio 12.5), e = 0.9
  c = october_kite()
  p = october_kite(polar = :profile)
  induced = 1 / (pi * 0.9 * 26.0^2 / 54.0)
  @test induced_drag_factor(26.0, 54.0) ≈ induced
  for c_l in (0.5, 1.0, 2.0)
    @test calc_c_d(c, c_l, coeffs = c.c_d_fun_coeffs) - calc_c_d(p, c_l, coeffs = p.c_d_fun_coeffs) ≈ induced * c_l^2
  end
  @test 2.0 / calc_c_d(p, 2.0, coeffs = p.c_d_fun_coeffs) ≈ 116.7 rtol = 0.001
  @test 2.0 / calc_c_d(c, 2.0, coeffs = c.c_d_fun_coeffs) ≈ 15.37 rtol = 0.001
  @test_throws ArgumentError october_kite(polar = :other)
end


@testset "The shaft carries Lambda at its twist" begin
  (r0, r1, l) = (6.0, 100.0, 250.0)
  lambda_of(t) = sin(t) * r0 / sqrt(2 * cos(t) * r0 * r1 + l^2 - r0^2 - r1^2) # eq. (14)
  theta = shaft_twist(r0, r1, l, 0.025)
  @test lambda_of(theta) ≈ 0.025 rtol = 1e-9
  @test 0 < theta < pi / 2
  # 6 m carries 0.0262 at most, 5.73 m 0.025
  @test isnan(shaft_twist(r0, r1, l, 0.027))
  @test !isnan(shaft_twist(r0, r1, l, 0.026))
  @test isnan(shaft_twist(5.7, r1, l, 0.025))
  @test !isnan(shaft_twist(5.75, r1, l, 0.025))
  @test_throws ArgumentError shaft_geometry(october_kite(), 0.03)
  @test shaft_geometry(october_kite(r0 = 0.0), 0.03).cartwheel == 1.0
end


@testset "Cartwheel factor" begin
  for b in (0.06, 0.2, 1 / 3)
    # untwisted, the tether runs straight out from r0 to r1
    @test cartwheel_factor(b, 0.0) ≈ (1 + b) * (1 + b^2) rtol = 1e-9
  end
  @test cartwheel_factor(0.0, 1.0) == 1.0
  # a small cartwheel: 1 + b cos(theta) + b^2 (1 + cos(theta)^2) / 2, to O(b^3)
  for theta in (0.0, 0.5, 1.2, pi / 2)
    b = 0.01
    @test cartwheel_factor(b, theta) - (1 + b * cos(theta) + b^2 * (1 + cos(theta)^2) / 2) ≈ 0.0 atol = 3e-6
  end
end


@testset "Tether drag against the reference model" begin
  # The arm-weighted drag over D_full = rho C_t d l v^2 / 2 from the reference model of the author's
  # tether-shape work (steady whirling tether, variable tension, centrifugal force, cross-flow drag,
  # no gravity, no wind), solved by shooting: tether-check/reference_check.jl of the revised paper,
  # with the ground end on the axis or on a cartwheel of 6 m at the twist of Lambda = 0.025.
  # October tether: 40 mm, 250 m, C_t 1.1, kite at 100 m; v and the tension along the tether at the
  # kite as below.
  cases = [(r0 = 0.0, v = 100.0, tension = 760e3, reference = 0.252508, rtol = 1e-4),
           (r0 = 0.0, v = 120.0, tension = 150e3, reference = 0.264512, rtol = 5e-4),
           (r0 = 6.0, v = 100.0, tension = 760e3, reference = 0.256873, rtol = 3e-3)]
  for k in cases
    c = october_kite(r0 = k.r0)
    shaft = shaft_geometry(c, 0.025)
    # the tension along the shaft that gives this tension along the tether
    shaft_tension = k.tension * shaft.height / c.l
    ratio = tether_c_d(c, shaft, k.v, k.v, shaft_tension) * c.s / (c.c_d_tether * c.d * c.l)
    @test ratio ≈ k.reference rtol = k.rtol
    # a straight tether from the axis, 1/4, is 1 to 5 % below
    @test 0.25 < k.reference < 0.27
  end
  # the shape factor is 1 for a massless tether with no drag slope, and the belly alone follows
  # the closed form (cos 3 Kl - 9 cos Kl + 8) / (3 Kl sin(Kl)^3) to fourth order
  @test tether_shape_factor(0.0, 0.0, 0.4) == 1.0
  kl = 0.3
  @test tether_shape_factor(0.0, kl, 0.0) ≈ (cos(3kl) - 9cos(kl) + 8) / (3kl * sin(kl)^3) rtol = 1e-5
end


@testset "The simulation uses the tether's drag" begin
  c = october_kite()
  (df, input) = solve_sector_df(c, 12.0, 0.0, 0.025, 0.0, shaft_tension = 500e3)
  shaft = shaft_geometry(c, 0.025)
  @test df.c_d ≈ [calc_c_d(c, cl, coeffs = c.c_d_fun_coeffs) for cl = df.c_l] .+ df.c_d_tether
  @test df.c_d_tether ≈ [tether_c_d(c, shaft, v_k, abs(v_a), 500e3) for (v_k, v_a) = zip(df.v_k, df.v_a)]
  # about the nominal C_t d l / (4 S), with the cartwheel and the shape a few per cent more
  nominal = c.c_d_tether * c.d * c.l / (4 * c.s)
  @test 1.0 < mean(df.c_d_tether .* df.v_a .^ 2 ./ df.v_k .^ 2) / nominal < 1.05
end


@testset "Induction from the thrust coefficient" begin
  for a in (0.0, 0.1, 0.3, 0.39)
    @test induction_from_thrust_coefficient(4 * a * (1 - a)) ≈ a
  end
  # Buhl's line joins momentum theory at C_T = 0.96, a = 0.4
  @test induction_from_thrust_coefficient(0.96 + 1e-9) ≈ 0.4 atol = 1e-6
  @test induction_from_thrust_coefficient(2.0) ≈ 1.0
end
