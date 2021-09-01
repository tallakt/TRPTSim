using Plots
using StackedPlot


function plot_from_args_helper(init, args, add_plot_fun)
  result = init

  for arg in args
    if typeof(arg) == Pair{String}{PowerCurve}
      result = add_plot_fun(result, arg.first, arg.second)
    elseif typeof(arg) == Tuple{String,PowerCurve}
      result = add_plot_fun(result, arg[1], arg[2])
    elseif typeof(arg) == PowerCurve
      result = add_plot_fun(result, "", arg)
    end
  end

  result
end


function plot_power_curves(args...)
  #pc::PowerCurve)
  init = plot(xlabel = "wind [m/s]", ylabel = "power of all kites [kW]")

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

  plot_from_args_helper(init, args, add_plot)
end


function plot_tension_curves(args...)
  #pc::PowerCurve)
  init = plot(xlabel = "wind m/s", ylabel = "tension ton")

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

  plot_from_args_helper(init, args, add_plot)
end


function plot_vmg_curves(args...)
  #pc::PowerCurve)
  init = plot(xlabel = "wind m/s", ylabel = "VMG")

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

  plot_from_args_helper(init, args, add_plot)
end


function heatmap_tension_moment_power(c::Configuration, wind::Number, psi::Number, shaft_tensions, mtrs, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = heuristic_flying_speed(c, wind, psi))

  function fun(t, mtr)
    (df,_)=TRPTSim.solve_sector_df(c, wind, psi, t, mtr, force_h, step_size = step_size, iterations = iterations, finish_threshold = finish_threshold, speed0 = speed0) 
    if maximum(df.c_l) > c.design_c_l
      NaN
    else
      mean(df.power) * c.n
    end
  end

  heatmap(mtrs
          , shaft_tensions
          , [fun(t, mtr) for t = shaft_tensions, mtr=mtrs]
         )
end


function plot_solution(df::DataFrame, solver_input::Dict)
  c = solver_input[:config]

  plotstacked(df.psi_p .|> rad2deg
              , trace([df.power signal_sum_of_kites(df.power, c.n)] .* 0.001 , lab="power [kW]", ylims=(0,Inf))
              , trace([df.moment df.sum_moments], lab=["moment [Nm]" "sum moments"], ylims=(0, Inf))
              , trace(df.v_k, lab="v_k [m/s]", ylims=(0,maximum(df.v_k) * 1.1))
              , trace(df.phi_p .|> rad2deg, lab="phi_p [deg]")
              , trace([df.shaft_tension signal_sum_of_kites(df.shaft_tension, c.n)] ./ 9.81 .* 0.001, lab=["shaft t [ton]" "sum"], ylims=(0,Inf))
              , trace([ones(size(df.c_l)) .* c.design_c_l df.c_l], lab=["" "C_L"], lc=[:lightgray 1], ylims=(0, c.design_c_l * 1.1))
              , trace([df.force_h df.force_v signal_sum_of_kites(df.force_h, c.n) signal_sum_of_kites(df.force_v, c.n)] ./ 9.81, lab=["force h [kg]" "force v" "sum force h" "sum force v"])
              , xlims=(0, 500)
             )
end
