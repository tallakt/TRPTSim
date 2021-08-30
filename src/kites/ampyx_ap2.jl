function ampyx_ap2()
  # From  Elena Malz, Jonas Koenemann, Sören Sieberling, Sébastien Gros (2019).
  # A reference model for airborne wind energy systems for optimization and
  # control. Renewable Energy.
  # http://awesco.eu/publication/malz-2019-a/malz-2019-a.pdf
  #
  # M/T factor 1-2
  config(n = 3
         , m = 36.8
         , s = 3.0    # b = 5.5, c = 0.55
         , d = 0.0030
         , l = 60.0
         , design_c_l = 1.4
         # I used the C_{X_0} and C_{Z_0} to arrive at this curve, but had to
         # add 0.1 to the drag coeff to not get negative drag at low C_L
         , c_d_fun_coeffs = [0.15479410788134904, -2.42091569028119,
                             19.374252384596552, -74.17189470394626,
                             155.68585189717493, -187.80138066307026,
                             129.8435082133324, -47.77113556476243,
                             7.248344782408041]
         , elev = deg2rad(30.0)
         , radius = 30.0
        )
end



