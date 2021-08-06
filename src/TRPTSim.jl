module TRPTSim

export config, config_to_dict, scale, scale_to_area, calc_c_d,
c_d_coeffs_with_tether, wind_vec, rotary_inertia, solve_sector,
solve_sector_df, optimal_tension_mtf, optimal_tension, shaft_section_c_d,
shaft_section_m_t_factor, shaft_section_compression, power_curve,
plot_power_curves, heatmap_tension_moment_power, heuristic_flying_speed,
heuristic_shaft_tension

include("configuration.jl")
include("tether.jl")
include("simulation.jl")
include("find_params.jl")
include("shaft.jl")
include("power_curve.jl")
include("plots.jl")
include("heuristics.jl")
end
