using Statistics

struct PowerCurveEntry
    wind::Float64
    power::Float64
    shaft_tension::Float64
end


struct PowerCurve
    curve::Array{PowerCurveEntry}
    config::Configuration
    psi::Float64
    moment_tension_factor::Float64
    force_h::Float64
end




function power_curve(c::Configuration, winds, psi::Number, mtr::Number, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = 100.0, shaft_tension_fun::Function = w -> optimal_tension(c, w, psi, mtr, force_h))
    make_entry = function(wind::Number) 
        shaft_tension = shaft_tension_fun(wind)
        res = solve_sector_df(c, wind, psi, mtr, force_h, iterations = iterations, step_size = step_size, finish_threshold = finish_threshold, speed0 = speed0, shaft_tension = shaft_tension)
        PowerCurveEntry(wind, get_avg_power(res), shaft_tension * c.n)
    end
    PowerCurve(map(make_entry, winds), c, psi, mtr, force_h)
end

