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
