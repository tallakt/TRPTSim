module TRPTSim

export config, config_to_dict, scale, scale_to_area, calc_c_d,
c_d_coeffs_with_tether, tether_mass, tether_strength, wind_vec, rotary_inertia,
solve_sector, solve_sector_df, optimal_tension_mtf, optimal_tension,
shaft_section_c_d, shaft_section_m_t_factor, shaft_section_compression,
power_curve, plot_power_curves, plot_tension_curves,
heatmap_tension_moment_power, plot_solution, heuristic_flying_speed,
heuristic_shaft_tension_per_kite, grid_optimize_1d, get_avg_power,
signal_sum_of_kites, october_kite

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
end
