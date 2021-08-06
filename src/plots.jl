using Plots


function plot_power_curves(args...)
  #pc::PowerCurve)
  result = plot()

  function add_plot(p, label, pc)
    plot!(p
          , [entry.wind for entry=pc.curve]
          , [entry.power * 0.001 for entry=pc.curve]
          , ylims = (0, Inf)
          , xlims = (0, Inf)
          , lab = label
          , legend = :topleft
         )
  end

  for arg in args
    if typeof(arg) == Pair{String}{PowerCurve}
      result = add_plot(result, arg.first, arg.second)
    elseif typeof(arg) == PowerCurve
      result = add_plot(result, "", arg)
    end
  end

  result
end


function heatmap_tension_moment_power(c::Configuration, wind::Number, psi::Number, shaft_tensions, m_t_factors, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = 100.0)

  function fun(t, m_t)
    (df,_)=TRPTSim.solve_sector_df(c, wind, psi, t, m_t, force_h, step_size = step_size, iterations = iterations, finish_threshold = finish_threshold, speed0 = speed0) 
    if maximum(df.c_l) > c.design_c_l
      NaN
    else
      mean(df.power)
    end
  end

  display(MIME"image/png"()
          , heatmap(m_t_factors
                    , shaft_tensions
                    , [fun(t, m_t) for t = shaft_tensions, m_t=m_t_factors]
                   )
         )
end
