"""
    m600_kite(; e = 0.85, c_l = 2.56, c_d = 0.2443, design_c_l = c_l,
              radius = 145.0, l = 340.0, r0 = 9.0, d = 0.004 * sqrt(3 * 250 / 21))

Three Makani M600 kites as built, the kite that flew, on a plant sized for the M600 by the
heuristics of the revised Pyramid paper: a ring of `radius`, the M600's loop radius in the bulk of
its flights, 145 m; a shaft of `l`, the shortest at which the lowest wing tip clears the cartwheel
by 17 m at the elevation of 30 degrees, 340 m; a cartwheel of `r0`, 9 m, kept small for the
tower; and a tether of diameter `d` for the M600's own design tension, 250 kN (Makani, Part I,
sec. 1.3.4), at the safety factor 3 (23.9 mm).

The kite is SN6 as built (Makani, The Energy Kite, Part I, "M600 Energy Kite Description",
Table 1): span 25.66 m, wing area 32.9 m^2, so an aspect ratio of 20.0, and a mass of 1730.8 kg.
It is a fly-gen kite, and it carries its rotors and pylons here, as mass and as drag.

The polar is C_D = C_D0 + C_L^2 / (pi e AR), with the span efficiency `e` (none is published) and
C_D0 such that C_D = `c_d` at C_L = `c_l`: by default Makani's as-built kite drag coefficient,
0.2443 at the lift coefficient of its best performance, 2.56 (Part I, "Airborne Wind Turbine
Performance", Table 1, "M600 As-Built"). `design_c_l`, the highest C_L the kites fly, is that C_L.
"""
function m600_kite(; e::Number = 0.85, c_l::Number = 2.56, c_d::Number = 0.2443, design_c_l::Number = c_l,
                   radius::Number = 145.0, l::Number = 340.0, r0::Number = 9.0,
                   d::Number = 0.004 * sqrt(3 * 250 / 21))
    span = 25.66
    area = 32.9
    k = induced_drag_factor(span, area, e = e)
    c_d0 = c_d - k * c_l^2
    c_d0 > 0 || throw(ArgumentError("C_D = $c_d at C_L = $c_l is less than the induced drag at e = $e"))
    config(n = 3
           , m = 1730.8
           , s = area
           , d = d
           , l = l
           , design_c_l = design_c_l
           , c_d_fun_coeffs = [c_d0, 0.0, k]
           , elev = deg2rad(30.0)
           , radius = radius
           , r0 = r0
          )
end


"""
    m600_ground_gen_kite(; mass_removed = 0.35, drag_removed = 0.0101, kwargs...)

Three M600 kites rebuilt for ground generation, on the plant of `m600_kite`: the M600's airframe
as built without its generation hardware (rotors, motor-generators, nacelles, pylons, power
electronics). Rough estimates, for want of a published breakdown:
- the mass less the share `mass_removed` of 1730.8 kg: Makani writes that "~1/3 to 1/2 the total
  mass of the kite will be the power system" (The Energy Kite, Part I, "Airborne Wind Turbine
  Performance", sec. 7.1), its MX2 design roll-up has 35 % (pylons, nacelles with powertrains,
  wing tip thrusters; Part II, "The Spreadsheet Kite Structure Sizer"), and the M600's structure
  list has pylons 175 kg and rotors with spinners 65 kg, 14 % (Part II, "M600 Airframe");
- the kite drag coefficient less `drag_removed`, on C_D0: a build-up of the pylons (the sections of
  Makani's ASWING model of the M600), their junctions, the nacelles and the idling rotors with
  Hoerner's component drag coefficients, 0.0101 (0.0068 to 0.0171), against 0.2443 as built.
`kwargs` go to `m600_kite`: the polar's anchor and span efficiency, the cartwheel.
"""
function m600_ground_gen_kite(; mass_removed::Number = 0.35, drag_removed::Number = 0.0101, kwargs...)
    c = m600_kite(; kwargs...)
    coeffs = copy(c.c_d_fun_coeffs)
    coeffs[1] -= drag_removed
    config(c, m = c.m * (1 - mass_removed), c_d_fun_coeffs = coeffs)
end
