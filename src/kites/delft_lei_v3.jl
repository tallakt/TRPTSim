function delft_lei_v3()
  # https://wes.copernicus.org/articles/4/1/2019/
  # M/T factor 2.0 is fine for this kite
  config(n = 3
         , m = 22.8
         , s = 19.75 # 25.0 flat, wingspan 8.3 m, chord 2.2 m
         , d = 0.0065
         , l = 60.0
         , design_c_l = 1.5
         # Just make a linear C_L/C_D for G_e = 4.0
         , c_d_fun_coeffs = [0.01, 1 / 4.0]
         , elev = deg2rad(30.0)
         , radius = 30.0 # wingspan 42.7, Æ 12.1
        )
end


