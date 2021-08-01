module TRPTSim

export config, scale, scale_to_area, calc_c_d, wind_vec, rotary_inertia, solve_sector

include("configuration.jl")
include("tether.jl")
include("simulation.jl")

end
