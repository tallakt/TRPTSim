
DENSITY_OF_UHWMPE_TETHER = 711.5

function tether_mass(diameter, length = 1.0)
    (diameter / 2)^2 * pi * DENSITY_OF_UHWMPE_TETHER * length
end

function tether_strength(diameter)
    diameter ^ 2 / 0.004^2 * 21_000 # using DSK90 4 mm as reference
end


# Tether drag
#
# The drag of a tether costs the kites power through the moment it makes about the axis of the
# shaft, and the kite pays for it as if it had that moment over its looping radius r1 as extra drag
# (the arm weighting). For a straight tether from the axis this is a quarter of
# D_full = 1/2 rho C_t d l v_k^2, the drag the tether would have if all of it moved at the kite's
# speed v_k. Two factors refine it:
# - the cartwheel: the tether runs straight from the cartwheel, radius r0, twisted theta behind the
#   kite; cartwheel_factor;
# - the shape: centrifugal force bends the tether out into a belly, and its drag bends it back
#   into a spiral that lags the kite; tether_shape_factor.
# The shape factor is the series of the steady whirling tether (variable tension, centrifugal
# force, cross-flow drag, no gravity, no wind, ground end on the axis) to second order, which the
# author's tether-shape work derived and checked against a numerical solution of that model.


"""
    shaft_twist(r0, r1, l, lambda)

The twist (radians) at which a shaft of one section, cartwheel radius `r0`, looping radius `r1` and
tether length `l`, carries the moment `lambda` (Lambda, moment per tension per looping radius): the
twist below that of the shaft's largest Lambda. NaN if the shaft can't carry `lambda`.
"""
function shaft_twist(r0::Number, r1::Number, l::Number, lambda::Number)
    lambda <= 0 && return 0.0
    r0 <= 0 && return NaN
    h_0 = sqrt(l^2 - (r1 - r0)^2)
    h_pi = sqrt(l^2 - (r1 + r0)^2)
    theta_max = acos((h_0 * h_pi - (l^2 - r0^2 - r1^2)) / (2 * r0 * r1))
    lambda_of(t) = sin(t) * r0 / sqrt(2 * cos(t) * r0 * r1 + l^2 - r0^2 - r1^2)
    lambda > lambda_of(theta_max) && return NaN
    (a, b) = (0.0, theta_max)
    for _ in 1:60
        m = (a + b) / 2
        lambda_of(m) < lambda ? (a = m) : (b = m)
    end
    (a + b) / 2
end


"""
    cartwheel_factor(b, theta)

The drag moment of a straight tether from a cartwheel of radius `b` r1, twisted `theta` behind the
kite, to the kite at r1, over that of a straight tether from the axis: 4 times the integral over
t = 0..1 of rho(t)^3, where rho(t) r1 is the tether's distance from the axis. Untwisted it is
(1 + b)(1 + b^2); for a small cartwheel about 1 + b cos(theta) + b^2 (1 + cos(theta)^2) / 2.
"""
function cartwheel_factor(b::Number, theta::Number)
    b == 0 && return 1.0
    f(t) = ((1 - t)^2 * b^2 + 2 * t * (1 - t) * b * cos(theta) + t^2)^1.5
    n = 64 # Simpson's rule; the integrand is smooth
    h = 1 / n
    4 * h / 3 * (f(0.0) + f(1.0) + 4 * sum(f((2 * i - 1) * h) for i in 1:(n ÷ 2)) +
                 2 * sum(f(2 * i * h) for i in 1:(n ÷ 2 - 1)))
end


"""
    tether_shape_factor(nu, kl, eps)

How much the shape of a whirling tether changes its arm-weighted drag, against a straight tether:
1 + order 1 + order 2 of the series in the drag slope at the kite `nu` = rho C_t d l v^2 / (8 T), the
belly number `kl` = K l with K^2 = mu omega^2 / T, and `eps` = r1 / l, with T the tension along the
tether at the kite. In sigma = nu / eps = nu l / r1, order 1 is
(K l)^2 / 6 - (K r1)^2 / 6 - 3 sigma^2 / 70 - 3 nu^2 / 7: the belly adds, the spiral pulls the tether
in and takes away. The series was checked against the numerical model for nu <= 0.2, sigma <= 2 and
K l <= 1; outside that range it is taken at the edge (states far from any optimum).
"""
function tether_shape_factor(nu::Number, kl::Number, eps::Number)
    nu = min(nu, 0.2, 2 * eps)
    kl = min(kl, 1.0)
    s2 = nu == 0 ? 0.0 : (nu / eps)^2
    k2 = kl^2
    e2 = eps^2
    order_1 = k2 / 6 * (1 - e2) - 3 * s2 / 70 - 3 * nu^2 / 7
    order_2 = e2^2 * (k2^2 / 18 + 65 * k2 * s2 / 252 + 4995 * s2^2 / 20384) +
              e2 * (-61 * k2^2 / 720 - 653 * k2 * s2 / 3780 + 3057 * s2^2 / 50960) +
              7 * k2^2 / 240 - 5 * k2 * s2 / 189 + 411 * s2^2 / 101920
    1 + order_1 + order_2
end


"""
    shaft_geometry(c, mtr)

The shaft of `c` carrying Lambda = `mtr`: its twist, height (along the axis) and cartwheel factor.
With `c.r0 = 0` the tethers come from the axis and the twist is not computed.
"""
function shaft_geometry(c::Configuration, mtr::Number)
    if c.r0 <= 0
        return (theta = 0.0, height = sqrt(c.l^2 - c.radius^2), cartwheel = 1.0)
    end
    theta = shaft_twist(c.r0, c.radius, c.l, mtr)
    isnan(theta) && throw(ArgumentError("a cartwheel of r0 = $(c.r0) m can't carry Lambda = $mtr"))
    height = sqrt(c.l^2 - c.r0^2 - c.radius^2 + 2 * c.r0 * c.radius * cos(theta))
    (theta = theta, height = height, cartwheel = cartwheel_factor(c.r0 / c.radius, theta))
end


"""
    tether_c_d(c, shaft, v_k, v_a, shaft_tension)

The drag of one tether as a drag coefficient of its kite, on the kite's area and its apparent speed
`v_a`, for the kite flying at `v_k` (= omega r1) with the tension `shaft_tension` along the shaft:
C_t d l / (4 S) (v_k / v_a)^2 times the cartwheel factor of the shaft (`shaft_geometry`) and the
shape factor.
"""
function tether_c_d(c::Configuration, shaft, v_k::Number, v_a::Number, shaft_tension::Number)
    tension = max(shaft_tension, 1.0) * c.l / shaft.height # along the tether
    nu = c.rho * v_k^2 * c.c_d_tether * c.d * c.l / (8 * tension)
    kl = v_k / c.radius * sqrt(tether_mass(c.d) / tension) * c.l
    shape = tether_shape_factor(nu, kl, c.radius / c.l)
    c.c_d_tether * c.d * c.l / (4 * c.s) * (v_k / v_a)^2 * shaft.cartwheel * shape
end
