using Roots
using Plots
using Rotations
using LinearAlgebra
using LaTeXStrings


function lambda_helper(r0, r1, l, theta)
  sin(theta) * r0 / sqrt(2 * cos(theta) * r0 * r1 + l^2 - r0^2 - r1^2 + 0im)
end


function lambda(r0, r1, l, theta)
  tmp = lambda_helper(r0, r1, l, theta)
  isreal(tmp) ? real(tmp) : NaN
end


function lambda_error(r0, r1, r2, l0, l1, theta0, theta1)
  lambda0 = lambda_helper(r0, r1, l0, theta0) * r1 / r2
  lambda1 = lambda_helper(r1, r2, l1, theta1)
  (real(lambda0 - lambda1))^2 + imag(lambda0)^2 + imag(lambda1)^2
end


function relax_twist(r0, r1, r2, l0, l1, theta)
  result = find_zero(theta0 -> lambda_error(r0, r1, r2, l0, l1, theta0, theta - theta0), theta * l0 / (l0 + l1))
end


function find_double_section_lambda(r0, r1, r2, l0, l1, theta)
  # like before, but check if tether force is pointing in or out by using
  # vector math
  #
  # - calc vectors
  # - calc tensions
  # - calc sum of forces
  # - dot product with outward vector
  # - bigger than 0 -> ok, else NaN
  # point are x0, x1, x2
  try
    theta0 = relax_twist(r0, r1, r2, l0, l1, theta)
    theta1 = theta - theta0
    x0 = RotZ(-theta0) * [r0; 0; 0]
    x1 = [r1; 0; 0]
    x2 = RotZ(theta1) * [r2; 0; 0]

    # assume T_N = 1.0
    tension0 = l0 / sqrt(2 * cos(theta0) * r0 * r1 + l0^2 - r0^2 - r1^2) 
    tension1 = l1 / sqrt(2 * cos(theta1) * r1 * r2 + l1^2 - r1^2 - r2^2) 

    force0 = normalize(x0 - x1) * tension0
    force1 = normalize(x2 - x1) * tension1

    radial_force = dot(force0 + force1, [1.0; 0.0; 0.0])

    if isreal(tension0) && isreal(tension1) && radial_force >= 0
      lambda(r1, r2, l1, theta - theta0)
    else
      # soft bridle acts compressively, or invalid geometry
      NaN
    end
  catch
    NaN
  end
end


function find_possible_buckling_theta(r0, r1, r2, l0, l1)
  try
    (true, find_zero(th -> lambda(r0, r2, l0 + l1) - find_double_section_lambda(r0, r1, r2, l0, l1, th), pi/2))
  catch
    (false, NaN)
  end
end


# find the point where the two section shaft buckes. This should be the same
# point where the single segment shaft without the bridle has the same lambda
# value
function max_lambda(r0, r1, r2, l0, l1; accuracy = 0.001)
  (buckling, theta) = find_possible_buckling_theta(r0, r1, r2, l0, l1)
  theta_max = buckling ? theta : deg2rad(180)
  res = if buckling
    lambda(r0, r2, l0 + l1, theta)
  else
    tmp = maximum(filter(x -> !isnan(x), vcat(-1, find_double_section_lambda.(r0, r1, r2, l0, l1, 0:accuracy:theta_max))))
  end
  res >= 0 ? res : NaN
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
