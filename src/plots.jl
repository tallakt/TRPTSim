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
