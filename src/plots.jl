using Plots
using StackedPlot


function plot_power_curves(args...)
  #pc::PowerCurve)
  result = plot(xlabel = "wind m/s", ylabel = "power kW")

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


function plot_tension_curves(args...)
  #pc::PowerCurve)
  result = plot(xlabel = "wind m/s", ylabel = "tension ton")

  function add_plot(p, label, pc)
    plot!(p
          , [entry.wind for entry=pc.curve]
          , [entry.shaft_tension ./ 9.81 * 0.001 for entry=pc.curve]
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


function plot_vmg_curves(args...)
  #pc::PowerCurve)
  result = plot(xlabel = "wind m/s", ylabel = "VMG")

  function add_plot(p, label, pc)
    plot!(p
          , [entry.wind for entry=pc.curve]
          , [entry.power / entry.shaft_tension for entry=pc.curve]
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


function heatmap_tension_moment_power(c::Configuration, wind::Number, psi::Number, shaft_tensions, m_t_factors, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = heuristic_flying_speed(c, wind, psi))

  function fun(t, m_t)
    (df,_)=TRPTSim.solve_sector_df(c, wind, psi, t, m_t, force_h, step_size = step_size, iterations = iterations, finish_threshold = finish_threshold, speed0 = speed0) 
    if maximum(df.c_l) > c.design_c_l
      NaN
    else
      mean(df.power) * c.n
    end
  end

  heatmap(m_t_factors
          , shaft_tensions
          , [fun(t, m_t) for t = shaft_tensions, m_t=m_t_factors]
         )
end


function plot_solution(df::DataFrame, solver_input::Dict)
  c = solver_input[:config]

  plotstacked(df.psi_p .|> rad2deg
              , trace(df.power*0.001, lab="power", ylims=(0,Inf))
              , trace([df.moment df.sum_moments], lab=["moment" "sum moments"])
              , trace(df.v_k, lab="v_k", ylims=(0,Inf))
              , trace(df.phi_p .|> rad2deg, lab="phi_p")
              , trace(df.shaft_tension, lab="sh tension", ylims=(0,Inf))
              , trace(df.c_l, lab="C_L", ylims=(0,c.design_c_l))
              , trace([df.force_h df.force_v], lab=["force h" "force v"])
              , xlims=(0, 450)
             )
end
