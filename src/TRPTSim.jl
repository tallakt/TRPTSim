module TRPTSim

export config, config_to_dict, scale, scale_to_area, calc_c_d,
c_d_coeffs_with_tether, wind_vec, rotary_inertia, solve_sector,
solve_sector_df, optimal_tension_mtf, shaft_section_c_d,
shaft_section_m_t_factor, shaft_section_compression, power_curve ,
plot_power_curves

include("configuration.jl")
include("tether.jl")
include("simulation.jl")
include("find_params.jl")
include("shaft.jl")
include("power_curve.jl")
include("plots.jl")

# TODO power curve
# TODO quick parameter selector, and generate function
end
