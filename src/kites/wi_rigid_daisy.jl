function wi_rigid_daisy()
  # From Oliver Tulloch's paper
  # https://forum.awesystems.info/uploads/default/original/1X/b929b7041ed1a494b8f82574eae6abd2e3e56a3d.pdf
  #
  # Rigid Rotor
  # The rigid rotor set up has three foam wings that have NACA 4412 airfoil
  # profiles. These wings have a span of 1m and a chord of 0.2m. Similar to the
  # soft rotor the wings are equally spaced around a 4mm carbon fibre ring that
  # has a diameter of 3.04m. The outer tips of the wings are at a radius of
  # 2.22m therefore 0.28m of the wings span are located inside the ring. The
  # rigid rotor is designed to be as similar as possible to the soft rotor.
  # This allows for a more reliable comparison between the performance of the
  # HQ Kites and the foam wings.
  #
  # An M/(T R) factor of 0.05 is appropriate for this rig
  config(n = 3
         , m = 0.420
         , s = 1.0 * 0.2
         , d = 0.0012
         , l = 6.7
         , design_c_l = 0.635
         # airfoiltools polar Re 200.000 adjusted for AR 5.0
         # http://www.airfoiltools.com/polar/details?polar=xf-naca4412-il-200000-n5
         , c_d_fun_coeffs = [0.014820822624367794, -0.018815686974794157, 0.038103591224841526, 0.07548823309303064, 0.14556784366635053, -0.2614280022987595, -0.07065396953612338, 0.20733992500969514, -0.05974806561647341]
         , elev = deg2rad(30.0)
         , radius = 3.04 + 0.5 * (2.22 + 0.28)
        )
end


