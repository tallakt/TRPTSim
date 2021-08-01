MASS_SCALING_EXPONENT = 2.4


struct Configuration
  n :: Int                # number of kites
  m :: Float64            # mass per kite
  s :: Float64            # wing area per kite
  d :: Float64            # tether diameter
  l :: Float64            # Nominal shaft length
  design_c_l :: Float64   # nominal max C_L to use
  c_d_fun_coeffs :: Array{Float64, 1}
  w :: Float64            # wind speed, always coming from "north"
  elev :: Float64         # elevation angle, in radians, theta
  radius :: Float64       # looping radius
  gravity :: Float64      # always the same
  rho :: Float64          # air density
end


Configuration() = Configuration(3, 5.0, 0.5, 50.0, 0.005, 1.2, [0, 0.1], 12.0, deg2rad(30), 5.0, 9.81, 1.225)



function config(c::Configuration = Configuration(); kwargs...)
  merge(c |> config_to_dict, Dict{Symbol, Any}(kwargs)) |> dict_to_config
end


function scale(c::Configuration, scale::Number)
  config(m = c.m * scale^MASS_SCALING_EXPONENT
         , s = c.s * scale^2
         , d = c.d * scale
         , l = c.l * scale
         , radius = c.radius * scale
        )
end


function scale_to_area(c::Configuration, s::Number)
  scale(c, sqrt(s / c.s))
end

function calc_c_d(c::Configuration, c_l::Number)
  sum(map(i -> c.c_d_fun_coeffs[i] * c_l^(i-1), 1:length(c.c_d_fun_coeffs)))
end


function rotary_inertia(c::Configuration)
  c.n_kites * c.m * c.radius
end

# private

function config_to_dict(c::Configuration)
  foldl((acc, k) -> (acc[k] = getfield(c, k); acc), fieldnames(Configuration), init = Dict{Symbol, Any}())
end


function dict_to_config(dict::Dict{Symbol, Any})
  bad_keys = setdiff(keys(dict), fieldnames(Configuration))
  if length(bad_keys) > 0
    throw("Bad Keys: $(join(bad_keys, ", "))")
  end
  params = map(k -> dict[k], fieldnames(Configuration))
  Configuration(params...)
end
