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

  power_scaling = [
             ("W", 0, 3.0),
             ("kW", 3000.0, 0.001),
             ("MW", 3_000_000.0, 0.000_001),
            ];
  max_power = plot_from_args_helper(0, args, (_, _, x) -> maximum(map(y -> y.power, x.curve))) |> maximum
  (power_unit, _, power_factor) = power_scaling[findlast(x -> x[2] < max_power, power_scaling)]

  init = plot(xlabel = "wind [m/s]", ylabel = "power of all kites [$(power_unit)]")

  function add_plot(p, label, pc)
    plot!(p
          , [entry.wind for entry=pc.curve]
          , [entry.power * power_factor for entry=pc.curve]
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
  tension_scaling = [
             ("N", 0, 3.0),
             ("kN", 3000.0, 0.001),
             ("MN", 3_000_000.0, 0.000_001),
            ];
  max_tension = plot_from_args_helper(0, args, (_, _, x) -> maximum(map(y -> y.shaft_tension, x.curve))) |> maximum
  (tension_unit, _, tension_factor) = tension_scaling[findlast(x -> x[2] < max_tension, tension_scaling)]

  init = plot(xlabel = "wind m/s", ylabel = "shaft tension [$(tension_unit)]")

  function add_plot(p, label, pc)
    plot!(p
          , [entry.wind for entry=pc.curve]
          , [entry.shaft_tension .* tension_factor for entry=pc.curve]
          , ylims = (0, 1.1 * max_tension * tension_factor)
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

  power_scaling = [
             ("W", 0, 1.0),
             ("kW", 1000.0, 0.001),
             ("MW", 1_000_000.0, 0.000_001),
            ];
  (power_unit, _, power_factor) = power_scaling[findlast(x -> x[2] < maximum(df.power), power_scaling)]

  moment_scaling = [
             ("Nm", 0, 1.0),
             ("kNm", 1000.0, 0.001),
             ("MNm", 1_000_000.0, 0.000_001),
            ];
  (moment_unit, _, moment_factor) = moment_scaling[findlast(x -> x[2] < maximum(df.moment), moment_scaling)]

  force_scaling = [
             ("N", 0, 1.0),
             ("kN", 1000.0, 0.001),
             ("MN", 1_000_000.0, 0.000_001),
            ];
  (force_unit, _, force_factor) = force_scaling[findlast(x -> x[2] < maximum(df.force_v), force_scaling)]

  tension_scaling = [
             ("N", 0, 1.0),
             ("kN", 1000.0, 0.001),
             ("MN", 1_000_000.0, 0.000_001),
            ];
  (tension_unit, _, tension_factor) = tension_scaling[findlast(x -> x[2] < maximum(df.shaft_tension), tension_scaling)]

  plotstacked(df.psi_p .|> rad2deg
              , trace([df.power signal_sum_of_kites(df.power, c.n)] .* power_factor , lab=["power [$(power_unit)]" "sum power"], ylims=(0, 1.1 * maximum(signal_sum_of_kites(df.power, c.n) .* power_factor)))
              , trace([df.moment df.sum_moments] .* moment_factor, lab=["moment [$(moment_unit)]" "sum moments"], ylims=(0, 1.1 * maximum(df.sum_moments .* moment_factor)))
              , trace(df.v_k, lab="v_k [m/s]", ylims=(0,maximum(df.v_k) * 1.1))
              , trace(df.phi_p .|> rad2deg, lab="phi_p [deg]")
              , trace([df.shaft_tension signal_sum_of_kites(df.shaft_tension, c.n)] .* tension_factor, lab=["shaft tension [$(tension_unit)]" "sum"], ylims=(0, 1.1 * maximum(signal_sum_of_kites(df.shaft_tension, c.n) .* tension_factor)))
              , trace([ones(size(df.c_l)) .* c.design_c_l df.c_l], lab=["" "C_L"], lc=[:gray 1], ylims=(0, c.design_c_l * 1.1))
              , trace([df.force_h signal_sum_of_kites(df.force_h, c.n) df.force_v signal_sum_of_kites(df.force_v, c.n)] .* force_factor, lab=["force horiz [$(force_unit)]" "sum force h" "force vert" "sum force v"])
              , xlims=(0, 500)
             )
end
