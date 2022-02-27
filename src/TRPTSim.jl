module TRPTSim

export config, config_to_dict, scale, scale_to_area, calc_c_d,
c_d_coeffs_with_tether, tether_mass, tether_strength, wind_vec, rotary_inertia,
solve_sector, solve_sector_df, optimal_tension_mtf, optimal_tension,
shaft_section_c_d, shaft_section_m_t_factor, shaft_section_mtr,
shaft_section_cone_angle, shaft_section_compression, power_curve,
plot_power_curves, plot_tension_curves, heatmap_tension_moment_power,
plot_solution, heuristic_flying_speed, heuristic_shaft_tension_per_kite,
grid_optimize_1d, get_avg_power, signal_sum_of_kites,
get_bridle_alpha_b_and_force, max_alpha_b, october_kite, eijkelhof,
wi_rigid_daisy, ampyx_ap2, delft_lei_v3

include("configuration.jl")
include("tether.jl")
include("simulation.jl")
include("find_params.jl")
include("shaft.jl")
include("power_curve.jl")
include("plots.jl")
include("heuristics.jl")
include("grid_optimize.jl")
include("utilities.jl")
include("kites/october.jl")
include("kites/eijkelhof.jl")
include("kites/wi_rigid_daisy.jl")
include("kites/ampyx_ap2.jl")
include("kites/delft_lei_v3.jl")
end
