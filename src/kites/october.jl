function october_kite()
  config(n = 3
         , m = 1852.0
         , s = 54.0
         , d = 0.04
         , l = 160.0
         , design_c_l = 2.0
         # the glide number seems very high
         , c_d_fun_coeffs = [0.00436647190519774, 0.016732036920022945, -0.012590655090473097, 0.004275035381646313, -0.0002828891163939005]
         , elev = deg2rad(30.0)
         , radius = 80.0 # wingspan 26 m
        )
end
