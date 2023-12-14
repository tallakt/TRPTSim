function heuristic_flying_speed(c::Configuration, wind::Number, psi::Number)
    # without any torque power extraction
    c.design_c_l / calc_c_d(c, c.design_c_l) * wind * cos(sqrt(c.elev^2 + psi^2))
end


function heuristic_shaft_tension_per_kite(c::Configuration, wind::Number, psi::Number)
    # tension per kite
    speed = heuristic_flying_speed(c, wind, psi)
    0.5 * c.rho * c.design_c_l * 0.85 * c.s * (speed^2 + wind^2)
end

