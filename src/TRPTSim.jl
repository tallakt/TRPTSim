module TRPTSim

export config, config_to_dict, scale, scale_to_area, calc_c_d,
c_d_coeffs_with_tether, wind_vec, rotary_inertia, solve_sector,
solve_sector_df, optimal_tension_moment

include("configuration.jl")
include("tether.jl")
include("simulation.jl")

end
