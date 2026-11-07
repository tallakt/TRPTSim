"""
    october_kite(; polar = :physical, r0 = 6.0)

Three kites loosely based on the Makani 600 kW "October kite": span 26 m, 54 m^2, on 250 m of
40 mm tether, looping at 100 m.

`polar` picks the kite's drag polar:
- `:physical` (the default): the October polar below as the profile drag, plus the induced drag
  C_L^2 / (pi e AR) of the wing, AR = span^2 / area = 12.5 and e = 0.9. Glide number 15.4 at
  C_L = 2.
- `:profile`: the October polar alone, glide number 117 at C_L = 2: it has no induced drag.

`r0` is the radius of the cartwheel. 6 m carries Lambda = 0.026 at most, 5 % more than the 0.025
of the paper's runs; 5.73 m is the least for 0.025.
"""
function october_kite(; polar::Symbol = :physical, r0::Number = 6.0)
    span = 26.0
    area = 54.0
    profile = [0.00436647190519774, 0.016732036920022945, -0.012590655090473097, 0.004275035381646313, -0.0002828891163939005]
    coeffs = copy(profile)
    if polar == :physical
        coeffs[3] += induced_drag_factor(span, area, e = 0.9)
    elseif polar != :profile
        throw(ArgumentError("polar is :physical or :profile, not $polar"))
    end
    config(n = 3
           , m = 1852.0
           , s = area
           , d = 0.04
           , l = 250.0
           , design_c_l = 2.0
           , c_d_fun_coeffs = coeffs
           , elev = deg2rad(30.0)
           , radius = 100.0 # wingspan 26 m
           , r0 = r0
          )
end
