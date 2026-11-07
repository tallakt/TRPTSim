using Plots
using LaTeXStrings

# The two-section ("trousers") shaft: tethers from the cartwheel at r0 to a bridle at r1 over l0,
# then on to the kites at r2 over l1. Lambda is relative to r2.
#
# Until 2026-11 the maximum was found by a search for the twist where the two-section shaft and a
# single section r0 -> r2 have the same Lambda ("buckling"). That search was a secant search from
# 90 degrees over all twists: the two are equal untwisted and repeat every turn, so once
# its call was fixed (0a456b2) it found twists of 180 to 3060 degrees for l0 between about 11 and
# 16 m and drew nonsense there, and it never found the real buckling below l0 = 8 m. It also
# measured the bridle's load with the whole tether tension along the horizontal run, where the
# horizontal part of the tension is what pulls on the joint. Below, the shaft is worked out at
# each twist from 0 to 180 degrees instead.


# Lambda of one section from radius r0 to r1, length l, twist theta, relative to r1; NaN where the
# section is too short for the twist
function lambda(r0, r1, l, theta)
  h2 = 2 * cos(theta) * r0 * r1 + l^2 - r0^2 - r1^2
  h2 > 0 ? sin(theta) * r0 / sqrt(h2) : NaN
end


section_height(r0, r1, l, theta) = sqrt(2 * cos(theta) * r0 * r1 + l^2 - r0^2 - r1^2)


# The twist of the lower section when the shaft is twisted theta in all: both sections carry the
# same moment, so the moment difference, which goes from <= 0 at theta0 = 0 to >= 0 at theta0 =
# theta, is bisected at its first change of sign. NaN if the sections cannot be built.
function split_twist(r0, r1, r2, l0, l1, theta; n = 400)
  g = theta0 -> lambda(r0, r1, l0, theta0) * r1 / r2 - lambda(r1, r2, l1, theta - theta0)
  xs = range(0.0, theta, length = n)
  gs = g.(xs)
  for i = 1:(n - 1)
    (a, b, ga, gb) = (xs[i], xs[i + 1], gs[i], gs[i + 1])
    (isnan(ga) || isnan(gb)) && continue
    ga == 0 && return a
    if ga * gb < 0
      for _ = 1:60
        m = (a + b) / 2
        gm = g(m)
        sign(gm) == sign(ga) ? (a = m; ga = gm) : (b = m)
      end
      return (a + b) / 2
    end
  end
  NaN
end


# The radial pull of the two tethers on the joint at r1, outwards positive, per unit of the axial
# tension T_N: the horizontal part of a tether's tension is T_N times its horizontal run over its
# height. The bridle can only pull the joint in, so it holds while this is >= 0.
function joint_pull(r0, r1, r2, l0, l1, theta0, theta1)
  (r0 * cos(theta0) - r1) / section_height(r0, r1, l0, theta0) +
    (r2 * cos(theta1) - r1) / section_height(r1, r2, l1, theta1)
end


# Lambda of the shaft at the total twist theta. While the tethers pull the joint outwards the bridle
# holds it at r1 and the shaft has two sections. Once they would pull it in, the bridle goes slack,
# the joint moves in and each tether runs straight from r0 to r2: the shaft has buckled into a
# single section. NaN where neither can be built.
function double_section_lambda(r0, r1, r2, l0, l1, theta)
  theta0 = split_twist(r0, r1, r2, l0, l1, theta)
  if !isnan(theta0) && joint_pull(r0, r1, r2, l0, l1, theta0, theta - theta0) >= 0
    return lambda(r1, r2, l1, theta - theta0)
  end
  # buckled: the straight tether must pass the joint's height inside r1
  s = l0 / (l0 + l1)
  rho = sqrt((1 - s)^2 * r0^2 + s^2 * r2^2 + 2 * s * (1 - s) * r0 * r2 * cos(theta))
  rho <= r1 ? lambda(r0, r2, l0 + l1, theta) : NaN
end


# the largest Lambda over twists from 0 to 180 degrees
function max_lambda(r0, r1, r2, l0, l1; accuracy = 0.001)
  values = filter(!isnan, double_section_lambda.(r0, r1, r2, l0, l1, accuracy:accuracy:pi))
  isempty(values) ? NaN : maximum(values)
end


function do_plot(; length = 20, r0 = 1, r1 = 1, r2 = 4, accuracy = 0.001)
  let ll = 0:0.1:length
    plot(ll
         , max_lambda.(r0, r1, r2, ll, length .- ll, accuracy = accuracy)
         , lab=""
         , xlims = (0, length)
         , ylims = (0, Inf)
         , xlabel = "section split point"
         , ylabel = L"max $\Lambda$ possible"
        )
  end
end
