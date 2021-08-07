function october_kite()
  config(n = 3
         , m = 1852.0
         , s = 54.0
         , d = 0.04
         , l = 160.0
         , design_c_l = 2.2
         , c_d_fun_coeffs = [ 0.03667225862042596, -0.04533919437945615, 0.02901837617776685, -0.009015832926789383, 0.001164392508052731 ]
         , elev = deg2rad(30.0)
         , radius = 80.0 # wingspan 26 m
        )
end
