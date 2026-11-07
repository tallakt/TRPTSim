using LinearAlgebra

# Rotor induction, an estimate: the kites slow the wind they fly in. Actuator-disc momentum theory
# for the thrust of all kites along the wind gives the axial induction a, and the kites are flown
# in the wind (1 - a) w. The disc is taken across the shaft and projected on the wind. Momentum
# theory assumes the thrust spread evenly over the disc, and says nothing of the wake of the kites
# or of the skew of the disc; the area chosen sets the answer: the disc out to the wing tips is
# the main estimate, the annulus that the wings sweep a pessimistic bound.


"""
    induction_from_thrust_coefficient(c_t)

The axial induction a of an actuator disc with the thrust coefficient `c_t`: momentum theory,
c_t = 4 a (1 - a), up to a = 0.4 (c_t = 0.96), and Buhl's empirical line above.
"""
function induction_from_thrust_coefficient(c_t::Number)
    c_t <= 0.96 && return 0.5 * (1 - sqrt(1 - max(c_t, 0.0)))
    # c_t = 8/9 + (4 F - 40/9) a + (50/9 - 4 F) a^2 with the tip loss F = 1
    (A, B, C) = (50 / 9 - 4, 4 - 40 / 9, 8 / 9 - c_t)
    (-B + sqrt(B^2 - 4 * A * C)) / (2 * A)
end


"""
    aerodynamic_thrust(df, c)

The time-averaged aerodynamic force of all kites along the wind, from a `solve_sector_df` result.
"""
function aerodynamic_thrust(df::DataFrame, c::Configuration)
    along_wind = [dot(l + d, SVector(1.0, 0.0, 0.0)) for (l, d) = zip(df.lift_vec, df.drag_vec)]
    # the rows are equally spaced in the phase, so in time they weigh dt = dpsi / omega
    sum(along_wind ./ df.omega) / sum(1 ./ df.omega) * c.n
end


"""
    induced_wind(c, wind, psi, mtr, force_h, disc_area)

The wind (1 - a) `wind` that the kites fly in, with the induction a of an actuator disc of area
`disc_area` across the shaft (projected on the wind here) that takes the kites' thrust. The
tension is optimized for each wind tried. Returns `(wind, a, c_t, shaft_tension)`, NaN where the
kites don't fly.
"""
function induced_wind(c::Configuration, wind::Number, psi::Number, mtr::Number, force_h::Number,
                      disc_area::Number; iterations::Integer = 40, tolerance::Number = 1e-4)
    area = disc_area * cos(c.elev) * cos(psi) # the disc's normal is the shaft
    a = 0.0
    for _ in 1:iterations
        w = (1 - a) * wind
        tension = optimal_tension(c, w, psi, mtr, force_h)
        isnan(tension) && return (wind = NaN, a = NaN, c_t = NaN, shaft_tension = NaN)
        (df, _) = solve_sector_df(c, w, psi, mtr, force_h, shaft_tension = tension)
        c_t = aerodynamic_thrust(df, c) / (0.5 * c.rho * area * wind^2)
        a_new = induction_from_thrust_coefficient(c_t)
        abs(a_new - a) < tolerance && return (wind = w, a = a, c_t = c_t, shaft_tension = tension)
        a = 0.5 * (a + a_new) # damped
    end
    (wind = NaN, a = NaN, c_t = NaN, shaft_tension = NaN)
end
