using Statistics

function optimal_tension_mtf(c::Configuration, wind::Number, psi::Number, force_h::Number; step_size::Number = deg2rad(2.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0 = heuristic_flying_speed(c, wind, psi))
    # TODO DOESNT WORK, I THINK MAYBE THE CURVE IS NOT SMOOTH ENOUGH
    small_tension = heuristic_shaft_tension_per_kite(c, wind, psi) * 0.1;
    m_t0 = 1.0

    minimizer = function(x::Vector)
        tension = x[1]
        m_t = x[2]
        moment = tension * m_t
        (df, _) = solve_sector_df(c, wind, psi, moment, force_h; step_size = step_size, iterations = iterations, speed0 = speed0, shaft_tension = tension)
        if size(df, 1) < 1 || any(df.c_l .>= c.design_c_l) || tension < 0 || m_t > 10.0
            # we dont want to go here
            NaN
        else
            -mean(df.power)
        end
    end
    opt = optimize(minimizer, [small_tension, m_t0],
                   Optim.Options(
                                 iterations = 300
                                 # g_tol = 1e-4
                                 # store_trace = true
                                 #, show_trace = true
                                 #, extended_trace = true
                                )
                  )
    (t, m_t) = (opt.minimizer[1], opt.minimizer[2])
    if opt.ls_success && t > 0.0
        (t, m_t)
    else
        (NaN, NaN)
    end
end


function optimal_tension(c::Configuration, wind::Number, psi::Number, mtr::Number, force_h::Number; step_size::Number = deg2rad(2.0), iterations::Integer = 50, finish_threshold::Number = 0.001, speed0 = heuristic_flying_speed(c, wind, psi))

    minimizer = function(tension::Float64)
        (df, _) = solve_sector_df(c, wind, psi, mtr, force_h; step_size = step_size, iterations = iterations, speed0 = speed0, shaft_tension = tension)
        if any(df.c_l .>= c.design_c_l) || tension < 0 || any(df.shaft_tension .> tether_strength(c.d) / c.safety_factor)
            # we dont want to go here
            NaN
        else
            if size(df, 1) > 0 
                if mtr > 0.0
                    -mean(df.power)
                else
                    # maximize tension if unloaded
                    -mean(df.shaft_tension)
                end
            else
                NaN
            end
        end
    end
    grid_optimize_1d(exp.(range(-4.0, 0.0, length = 30)) .* heuristic_shaft_tension_per_kite(c, wind, psi) * 3.0, minimizer)
end




