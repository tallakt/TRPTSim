function optimal_tension_mtf(c::Configuration, wind::Number, psi::Number, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0 = heuristic_flying_speed(c, wind, psi))
  small_tension = wind^2 * cos(psi)^2 * cos(c.elev)^2 * c.s;
  small_moment = small_tension / 4.0;
  
  minimizer = function(x::Vector)
    tension = x[1]
    moment = x[2]
    act_moment = clamp(moment, 0, 10 * tension) # dont want to go wild with moment, 1.5x is normally correct
    (df, _) = solve_sector_df(c, wind, psi, tension, moment, force_h; step_size = step_size, iterations = iterations, speed0 = speed0)
    if any(df.c_l .>= c.design_c_l) || tension < 0 || moment > tension * 20
      # we dont want to go here
      NaN
    else
      -sum(clamp.(df.power, 0, Inf))
    end
  end
  opt = optimize(minimizer, [small_tension, small_moment])
  (t, m) = (opt.minimizer[1], opt.minimizer[2])
  if opt.ls_success && t > 0.0
    (t, m)
  else
    (NaN, NaN)
  end
end


function optimal_tension(c::Configuration, wind::Number, psi::Number, m_t_factor::Number, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0 = heuristic_flying_speed(c, wind, psi))
  small_tension = wind^2 * cos(psi)^2 * cos(c.elev)^2 * c.s;
  
  minimizer = function(x::Vector)
    tension = x[1]
    moment = tension * m_t_factor
    (df, _) = solve_sector_df(c, wind, psi, tension, moment, force_h; step_size = step_size, iterations = iterations, speed0 = speed0)
    if any(df.c_l .>= c.design_c_l) || tension < 0
      # we dont want to go here
      NaN
    else
      -sum(clamp.(df.power, 0, Inf)) - (m_t_factor > 0 ? 0.0 : tension)
    end
  end
  opt = optimize(minimizer, [small_tension], BFGS())
  t = opt.minimizer[1]
  if opt.ls_success && t > 0.0
    t
  else
    NaN
  end
end


