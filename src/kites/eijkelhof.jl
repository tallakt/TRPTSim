function eijkelhof()
    # 1.2 MW design
    # https://repository.tudelft.nl/islandora/object/uuid%3Ad65cbee3-1546-4ce4-bfaf-bd3cb34923e3?collection=research
    parasitic_drag_coeff = 0.005
    config(n = 3
           , m = 6885.0
           , s = 150.3
           , d = 0.04
           , l = 350.0
           , design_c_l = 1.5
           # the glide number seems very high
           , c_d_fun_coeffs = [parasitic_drag_coeff, 1 / 18.0]
           , elev = deg2rad(30.0)
           , radius = 150.0 # wingspan 42.7, Æ 12.1
          )
end

