MASS_SCALING_EXPONENT = 2.4


struct Configuration
    n :: Int                # number of kites
    m :: Float64            # mass per kite
    s :: Float64            # wing area per kite
    d :: Float64            # tether diameter
    l :: Float64            # tether length, cartwheel to kite
    design_c_l :: Float64   # nominal max C_L to use
    c_d_fun_coeffs :: Array{Float64, 1}
    elev :: Float64         # elevation angle, in radians, theta
    radius :: Float64       # looping radius
    gravity :: Float64      # always the same
    rho :: Float64          # air density
    c_d_tether::Float64     # drag coefficient of tether
    safety_factor::Float64  # for the tether
    r0::Float64             # radius of the cartwheel at the ground; 0 puts the tethers on the axis
end


Configuration() = Configuration(3, 5.0, 0.5, 0.005, 50.0, 1.2, [0, 0.1], deg2rad(30), 5.0, 9.81, 1.225, 1.1, 3.0, 0.0)



function config(c::Configuration = Configuration(); kwargs...)
    merge(c |> config_to_dict, Dict{Symbol, Any}(kwargs)) |> dict_to_config
end


function scale(c::Configuration, scale::Number)
    config(c, m = c.m * scale^MASS_SCALING_EXPONENT
           , s = c.s * scale^2
           , d = c.d * scale
           , l = c.l * scale
           , radius = c.radius * scale
           , r0 = c.r0 * scale
          )
end


function scale_to_area(c::Configuration, s::Number)
    scale(c, sqrt(s / c.s))
end

function calc_c_d(c::Configuration, c_l::Number; coeffs = c_d_coeffs_with_tether(c))
    sum(map(i -> coeffs[i] * c_l^(i-1), 1:length(coeffs)))
end

function c_d_coeffs_with_tether(c::Configuration)
    # The kite's polar with the drag of a straight tether from the axis added, as a drag coefficient
    # of the kite: a piece of tether at radius r moves at v_k r / r1 and costs the kite its drag
    # times r / r1, in all 1/8 rho C_t d l v_k^2, which is C_t d l / (4 S) of the kite's
    # 1/2 rho S v_k^2. It does not depend on the kite's C_L, so it is the constant term of the polar
    # (tmp[1]); it was the C_L term (tmp[2]), which doubled it at C_L = 2. This nominal value is
    # used for start values only: the simulation uses tether_c_d (tether.jl), which adds the
    # cartwheel and the shape of the tether.
    tmp = zeros(size(c.c_d_fun_coeffs))
    tmp[1] = c.c_d_tether * c.d * c.l / 4.0 / c.s
    tmp + c.c_d_fun_coeffs
end


# The induced drag of a wing, C_D,i = k C_L^2, has k = 1 / (pi e AR), with the aspect ratio
# AR = span^2 / area and the span efficiency e.
induced_drag_factor(span::Number, area::Number; e::Number = 0.9) = 1 / (pi * e * span^2 / area)


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
