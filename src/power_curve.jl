using Statistics

struct PowerCurveEntry
    wind::Float64           # the free wind
    power::Float64
    shaft_tension::Float64
    induction::Float64      # a: the kites fly in (1 - a) wind; 0 without induction
end


struct PowerCurve
    curve::Array{PowerCurveEntry}
    config::Configuration
    psi::Float64
    moment_tension_factor::Float64
    force_h::Float64
end




# With a finite `disc_area` the kites fly in the wind that rotor induction leaves them
# (induced_wind, induction.jl); with the default they fly in the free wind.
function power_curve(c::Configuration, winds, psi::Number, mtr::Number, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = 100.0, shaft_tension_fun::Function = w -> optimal_tension(c, w, psi, mtr, force_h), disc_area::Number = Inf)
    make_entry = function(wind::Number)
        (w, a) = isfinite(disc_area) ? induced_wind(c, wind, psi, mtr, force_h, disc_area)[(:wind, :a)] : (wind, 0.0)
        isnan(w) && return PowerCurveEntry(wind, NaN, NaN, NaN)
        shaft_tension = shaft_tension_fun(w)
        res = solve_sector_df(c, w, psi, mtr, force_h, iterations = iterations, step_size = step_size, finish_threshold = finish_threshold, speed0 = speed0, shaft_tension = shaft_tension)
        PowerCurveEntry(wind, get_avg_power(res), shaft_tension * c.n, a)
    end
    PowerCurve(map(make_entry, winds), c, psi, mtr, force_h)
end

