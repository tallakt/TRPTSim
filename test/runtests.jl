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
             , m = 0.01                    # very small mass
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
             , m = 0.01                    # very small mass
             , radius = 15.0

           )
  w = 12.0
  psi = 0.0
  m_t_factor = 0.1
  force_h = 0.0
  tension_factor = 0.5
  # want to fly with half possible tension
  tension = tension_factor * (0.5 * c.rho * c.design_c_l * c.s * glide ^ 2 * w ^ 2)
  (df, solve_input) = solve_sector_df(c, w, psi, m_t_factor, force_h, shaft_tension = tension)

  @test mean(df.moment) ≈ tension .* m_t_factor                     atol = 50.0
  @test mean(df.power) ≈ tension .* m_t_factor .* mean(df.omega)    atol = 100.0
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
             , m = 0.01                    # very small mass
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
  expected_c_d = 0.4

  @test calc_c_d(c, c.design_c_l) ≈ expected_c_d                      atol = 0.02
  @test mean(df.v_a) ≈ norm([1, c.design_c_l / expected_c_d]) * w     atol = 2.0
end



@testset "Test elevation and phi" begin
  w = 12.0
  m_t_factor = 1.0
  force_h = 0.0
  c = config(eijkelhof_kite(), gravity = 0.0)
  elev_psi = [(0.0, 0.0), (30.0, 0.0), (60.0, 0.0), (0.0, 30.0), (0.0, 60.0), (30.0, 30.0)]
  powers = [solve_sector_df(config(c, elev = deg2rad(elv)), w, deg2rad(psi), m_t_factor, force_h) |> get_avg_power for (elv, psi) = elev_psi]

  @test powers[2] ≈ 1_400_000.0 rtol = 0.05 # not sure about the power values, at least they should be reducing with angle
  @test powers[3] ≈ 430_000.0   rtol = 0.10
  @test powers[4] ≈ powers[2]   rtol = 0.01
  @test powers[5] ≈ powers[3]   rtol = 0.01
  @test powers[6] ≈ 1_200_000   rtol = 0.05
end


@testset "Test scaling and effect of mass" begin
  w = 9.0
  m_t_factor = 1.0
  force_h = 0.0
  psi = 0.0
  c1 = eijkelhof_kite()
  c2 = config(scale_to_area(c1, c1.s / 2), m = c1.m / 2)
  c3 = scale_to_area(c1, c1.s / 2)
  power1 = solve_sector_df(c1, w, psi, m_t_factor, force_h) |> get_avg_power
  power2 = solve_sector_df(c2, w, psi, m_t_factor, force_h) |> get_avg_power
  power3 = solve_sector_df(c3, w, psi, m_t_factor, force_h) |> get_avg_power

  @test c3.m < c2.m # mass should scale more than x^2
  @test_broken power1 ≈ 2 * power2     rtol = 0.001
  @test power3 ≈ power2         rtol = 0.100
  @test power3 < power2 
end


@testset "Test horizontal and vertical force generation" begin
  w = 9.0
  m_t_factor = 1.0
  psi = 0.0
  tension = 300_000.0
  c = eijkelhof_kite()
  (sol1, _) = solve_sector_df(c, w, psi, m_t_factor, 0.0, shaft_tension = tension)
  (sol2, _) = solve_sector_df(c, w, psi, m_t_factor, 1000.0, shaft_tension = tension)
  (sol3, _) = solve_sector_df(c, w, psi, m_t_factor, 2000.0, shaft_tension = tension)
  (sol4, _) = solve_sector_df(c, w, psi, m_t_factor, -2000.0, shaft_tension = tension)

  f_v2 = signal_sum_of_kites(sol2.force_v, c.n) 
  @test minimum(f_v2) > c.m * c.gravity * 0.3 * c.n
  @test maximum(f_v2) < c.m * c.gravity * 1.8 * c.n
  @test_broken mean(f_v2) ≈ c.m * c.gravity * c.n    rtol = 0.05
  @test mean(sol1.force_h) ≈ 0.0      atol = 600.0
  @test mean(sol2.force_h) ≈ 1000.0   atol = 700.0
  @test mean(sol3.force_h) ≈ 2000.0   atol = 700.0
  @test mean(sol4.force_h) ≈ -mean(sol3.force_h) atol = 1500.0
end


