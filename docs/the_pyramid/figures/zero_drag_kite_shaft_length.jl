using Plots
using LaTeXStrings
using Roots


println("plot_for_3_mm_tether()")
function plot_for_3_mm_tether(; nondimensional = false, scale = 1.0, c_d_t = 1.1, d = 0.003, safety_factor = 4.0, nominal_wind = 9.0, elevation = 35.0, target_glide_ratio = 8.0)
  rho = 1.225;
  k_r = 3.0;
  r_0 = 8;
  r_1 = r_0 * k_r;
  glide_ratio_kite0 = 8.0;
  glide_ratio_kite1 = 12.0;
  glide_ratio_kite2 = 20.0;


  function belly_drag_shaft(rho, c_d_t, r_0, r_1, mu, d, l, omega, tension)
    -rho * c_d_t * (((r_0 ^ 2 + r_1 ^ 2) * cos(omega * sqrt(mu) * l * tension ^ (-1//2)) - 2 * r_0 * r_1) * sqrt(tension) * sin(omega * sqrt(mu) * l * tension ^ (-1//2)) + 2 * (r_0 * r_1 * cos(omega * sqrt(mu) * l * tension ^ (-1//2)) - r_0 ^ 2 / 2 - r_1 ^ 2 / 2) * l * sqrt(mu) * omega) * mu ^ (-1//2) * d * omega / sin(omega * sqrt(mu) * l * tension ^ (-1//2)) ^ 2 / 4
  end


  function strength_of_uhmwpe_tether(diameter)
    diameter ^ 2 / 0.004^2 * 21_000 # using DSK90 4 mm as reference 
  end


  function mass_per_meter_uhmwpe(d)
    DENSITY_OF_UHWMPE_TETHER = 711.4999626069974
    (d / 2)^2 * pi * DENSITY_OF_UHWMPE_TETHER
  end


  function calc_glide_ratio(rho, c_d_t, r_0, r_1, d, l, safety_factor, glide_ratio_kite)
    mu = mass_per_meter_uhmwpe(d)
    lift = strength_of_uhmwpe_tether(d) / safety_factor
    function objective_fun(glide_ratio)
      flying_speed = nominal_wind * glide_ratio * cos(deg2rad(elevation)) * 2 / 3 # 2/3 due to induction factor optimal
      omega = flying_speed / r_1
      tether_drag = belly_drag_shaft(rho, c_d_t, r_0, r_1, mu, d, l, omega, lift)
      kite_drag = lift / glide_ratio_kite
      glide_ratio - lift / (tether_drag + kite_drag)
    end
    tmp = find_zero(objective_fun, 15.0)
    tmp < 3 ? NaN : tmp
  end


  function plot_nondimensional()
    formatter_fun = x -> "$(Int(round(x / 1000)))k"
    let lengths = 60:300
      lift = strength_of_uhmwpe_tether(d * scale) / safety_factor
      glide_ratios = calc_glide_ratio.(rho, c_d_t, r_0 * scale, r_1 * scale, d * scale, lengths .* scale, safety_factor, Inf)
      glide_ratios_k0 = calc_glide_ratio.(rho, c_d_t, r_0 * scale, r_1 * scale, d * scale, lengths .* scale, safety_factor, glide_ratio_kite0)
      glide_ratios_k1 = calc_glide_ratio.(rho, c_d_t, r_0 * scale, r_1 * scale, d * scale, lengths .* scale, safety_factor, glide_ratio_kite1)
      glide_ratios_k2 = calc_glide_ratio.(rho, c_d_t, r_0 * scale, r_1 * scale, d * scale, lengths .* scale, safety_factor, glide_ratio_kite2)
      lambdas = 1 ./ 2 ./ glide_ratios 
      lambdas_k0 = 1 ./ 2 ./ glide_ratios_k0 
      lambdas_k1 = 1 ./ 2 ./ glide_ratios_k1 
      lambdas_k2 = 1 ./ 2 ./ glide_ratios_k2 
      lambda_max = r_0 ./ sqrt.(lengths.^2 .- r_0.^2 .- r_1.^2)
      plot(lengths ./ d
           , [lambdas lambdas_k0 lambdas_k1 lambdas_k2]
           , ylabel = L"\Lambda"
           , xlabel = L"tether length $\varrho = \frac{l}{d}$"
           , title = L"non-dimensional $\Lambda$ vs tether length, dimensioned for $w$ = %$(nominal_wind) m/s"
           , lab = ["zero kite drag" "tether + kite glide r $(glide_ratio_kite0)" "tether + kite glide r $(glide_ratio_kite1)" "tether + kite glide r $(glide_ratio_kite2)"]
           , legend = false
           , ylims = (0, 0.10)
           , xlims = (0, Inf)
           , size = (1000, 700)
           , xformatter = formatter_fun
           , lc = [:black :green :brown2 :royalblue]
          )
      plot!([0; maximum(lengths)] ./ d, [1; 1] ./ 2 ./ target_glide_ratio, lc = :gray, linestyle = :dash, lab = "glide ratio $(target_glide_ratio):1")
      plot!(lengths ./ d, lambda_max, lc = :red, linestyle = :dash, lab = L"maximum $\Lambda$")
    end
  end


  function plot_3mm()
    let lengths = 60:300
      lift = strength_of_uhmwpe_tether(d) / safety_factor
      glide_ratios = calc_glide_ratio.(rho, c_d_t, r_0, r_1, d, lengths, safety_factor, Inf)
      glide_ratios_k0 = calc_glide_ratio.(rho, c_d_t, r_0, r_1, d, lengths, safety_factor, glide_ratio_kite0)
      glide_ratios_k1 = calc_glide_ratio.(rho, c_d_t, r_0, r_1, d, lengths, safety_factor, glide_ratio_kite1)
      glide_ratios_k2 = calc_glide_ratio.(rho, c_d_t, r_0, r_1, d, lengths, safety_factor, glide_ratio_kite2)
      lambdas = 1 ./ 2 ./ glide_ratios 
      lambdas_k0 = 1 ./ 2 ./ glide_ratios_k0 
      lambdas_k1 = 1 ./ 2 ./ glide_ratios_k1 
      lambdas_k2 = 1 ./ 2 ./ glide_ratios_k2 
      lambda_max = r_0 ./ sqrt.(lengths.^2 .- r_0.^2 .- r_1.^2)
      plot(lengths
           , [lambdas lambdas_k0 lambdas_k1 lambdas_k2]
           , ylabel = L"\Lambda"
           , xlabel = L"tether length $l$"
           , title = L"$d$ = %$(Int(d * 1000)) mm, $r_0$ = %$(r_0) m and $r_1$ = %$(r_1) m, $w$ = %$(nominal_wind) m/s"
           , lab = ["zero kite drag" "tether + kite glide r $(glide_ratio_kite0)" "tether + kite glide r $(glide_ratio_kite1)" "tether + kite glide r $(glide_ratio_kite2)"]
           #, legend = false
           , ylims = (0, 0.10)
           , xlims = (0, Inf)
           , size = (1000, 700)
           , lc = [:black :green :brown2 :royalblue]
          )
      plot!([0; maximum(lengths)], [1; 1] ./ 2 ./ target_glide_ratio, lc = :gray, linestyle = :dash, lab = "glide ratio $(target_glide_ratio):1")
      plot!(lengths, lambda_max, lc = :red, linestyle = :dash, lab = L"maximum $\Lambda$")
    end
  end


  if (nondimensional) 
    plot_nondimensional()
  else
    plot_3mm()
  end
end


nothing
