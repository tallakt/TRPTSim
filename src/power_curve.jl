using Statistics

struct PowerCurveEntry
  wind::Float64
  power::Float64
  shaft_tension::Float64
  moment_tension_factor::Float64
end


struct PowerCurve
  curve::Array{PowerCurveEntry}
  config::Configuration
  psi::Float64
  force_h::Float64
end




function power_curve(c::Configuration, winds, psi::Number, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = 100.0, shaft_tension_moment_fun::Function = w -> optimal_tension_mtf(c, w, psi, force_h))
  make_entry = function(wind::Number) 
    (shaft_tension, m_t_factor) = shaft_tension_moment_fun(wind)
    (df, _) = solve_sector_df(c, wind, psi, shaft_tension, m_t_factor, force_h, iterations = iterations, step_size = step_size, finish_threshold = finish_threshold, speed0 = speed0)
    PowerCurveEntry(wind, mean(df.power), shaft_tension, m_t_factor)
  end
  PowerCurve(map(make_entry, winds), c, psi, force_h)
end


function tension_curve(c::Configuration, winds, psi::Number, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = 100.0, shaft_tension_moment_fun::Function = w -> optimal_tension_mtf(c, w, psi, force_h))
  make_entry = function(wind::Number) 
    (shaft_tension, m_t_factor) = shaft_tension_moment_fun(wind)
    (df, _) = solve_sector_df(c, wind, psi, shaft_tension, 0.0, force_h, iterations = iterations, step_size = step_size, finish_threshold = finish_threshold, speed0 = speed0)
    PowerCurveEntry(wind, mean(df.power), shaft_tension, m_t_factor)
  end
  PowerCurve(map(make_entry, winds), c, psi, force_h)
end
