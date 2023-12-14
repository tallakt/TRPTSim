using Plots
using StackedPlot
using LaTeXStrings


function plot_from_args_helper(init, args, add_plot_fun)
    result = init

    for arg in args
        if typeof(arg) == Pair{String}{PowerCurve} # eg using func("PC" => pc, ...)
            result = add_plot_fun(result, arg.first, arg.second)
        elseif typeof(arg) == Tuple{String,PowerCurve} # eg using func([("PC", pc), ("PC2", pc2)]...)
            result = add_plot_fun(result, arg[1], arg[2])
        elseif typeof(arg) == PowerCurve # eg using func(pc1, pc2, ...)
            result = add_plot_fun(result, "", arg)
        end
    end

    result
end


function plot_power_curves(args...; hawt = nothing)
    power_scaling = [
                     ("W", 0, 1.0),
                     ("kW", 3000.0, 0.001),
                     ("MW", 3_000_000.0, 0.000_001),
                    ];
    max_power = plot_from_args_helper(1.0, args, (acc, _, x) -> max(acc, maximum(map(y -> isnan(y.power) ? 0 : y.power, x.curve))))
    (power_unit, _, power_factor) = power_scaling[findlast(x -> x[2] < max_power, power_scaling)]
    println((max_power, power_unit, power_factor)) 

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

    # Prot HAWT as reference if they are supplied
    if !isnothing(hawt)
        for (n, g) = enumerate(groupby(hawt, [:model]))
            hawt_max_wind = 20.0
            # we scale to max avg power to compare power curves
            ii = g.wind .<= hawt_max_wind
            pp = [x > 0 ? x : NaN for x = g.power[ii]] .* max_power ./ maximum(g.power) * power_factor
            init = plot!(init, g.wind[ii], pp, lc = :lightsalmon, lab = "")
            j = findfirst(pp .>= (max_power * power_factor / 6 * (2 + n)))
            init = annotate!(init, g.wind[ii][j], pp[j], text(first(g.name), 8, valign = :bottom , halign = :right, color = :lightsalmon)) 
        end
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


function heatmap_tension_moment_power(c::Configuration, wind::Number, azimuth::Number, shaft_tensions, lambdas, force_h::Number; step_size::Number = deg2rad(1.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0::Number = heuristic_flying_speed(c, wind, azimuth))

    function fun(t, lambda)
        (df,_)=TRPTSim.solve_sector_df(c, wind, azimuth, t, lambda, force_h, step_size = step_size, iterations = iterations, finish_threshold = finish_threshold, speed0 = speed0) 
        if maximum(df.c_l) > c.design_c_l
            NaN
        else
            mean(df.power) * c.multiplicity
        end
    end

    heatmap(lambdas
            , shaft_tensions
            , [fun(t, lambda) for t = shaft_tensions, lambda=lambdas]
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
    (force_unit, _, force_factor) = force_scaling[findlast(x -> x[2] < maximum(abs.(df.force_v)), force_scaling)]

    tension_scaling = [
                       ("N", 0, 1.0),
                       ("kN", 1000.0, 0.001),
                       ("MN", 1_000_000.0, 0.000_001),
                      ];
    (tension_unit, _, tension_factor) = tension_scaling[findlast(x -> x[2] < maximum(abs.(df.shaft_tension)), tension_scaling)]

    plotstacked(df.psi_p .|> rad2deg
                , trace([df.power signal_sum_of_kites(df.power, c.multiplicity)] .* power_factor , lab=["power [$(power_unit)]" "sum power"], ylims=(0, 1.1 * maximum(signal_sum_of_kites(df.power, c.multiplicity) .* power_factor)))
                , trace([df.moment df.sum_moments] .* moment_factor, lab=["moment [$(moment_unit)]" "sum moments"], ylims=(0, 1.1 * maximum(df.sum_moments .* moment_factor)))
                , trace(df.v_k, lab="v_k [m/s]", ylims=(0,maximum(df.v_k) * 1.1))
                , trace(df.phi_p .|> rad2deg, lab="phi_p [deg]")
                , trace([df.shaft_tension signal_sum_of_kites(df.shaft_tension, c.multiplicity)] .* tension_factor, lab=["shaft tension [$(tension_unit)]" "sum"], ylims=(0, 1.1 * maximum(signal_sum_of_kites(df.shaft_tension, c.multiplicity) .* tension_factor)))
                , trace([ones(size(df.c_l)) .* c.design_c_l df.c_l], lab=["" "C_L"], lc=[:gray 1], ylims=(0, c.design_c_l * 1.1))
                , trace([df.force_h signal_sum_of_kites(df.force_h, c.multiplicity) df.force_v signal_sum_of_kites(df.force_v, c.multiplicity)] .* force_factor, lab=["force horiz [$(force_unit)]" "sum force h" "force vert" "sum force v"])
                , xlims=(0, 500)
               )
end


function drag_polar_plot(c::Configuration)
    cls = LinRange(0.0, c.design_c_l, 100)
    design_c_d = calc_c_d(c, c.design_c_l, coeffs = c.c_d_fun_coeffs)
    p0 = plot(xlims = (0, design_c_d * 1.1), ylims = (0, c.design_c_l * 1.1), xlabel = L"C_D", ylabel = L"C_L", legend = false)
    p1 = annotate!(p0, design_c_d * 0.1, c.design_c_l * 0.9, Plots.text(L"Design $C_L$: %$(c.design_c_l)", 9, :left))
    p2 = annotate!(p1, design_c_d * 0.1, c.design_c_l * 0.8, Plots.text("Glide ratio: $(round(c.design_c_l / design_c_d, digits = 1))", 9, :left))
    p3 = plot!(p2, [0, design_c_d * 2], [0, c.design_c_l * 2], lc = :gray, ls = :dash)
    p4 = plot!(p3, [calc_c_d(c, cl, coeffs = c.c_d_fun_coeffs) for cl = cls], cls, lc = :black)
    p5 = scatter!(p4, [design_c_d], [c.design_c_l], mc = :black)
end
