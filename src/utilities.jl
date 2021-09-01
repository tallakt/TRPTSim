using Statistics

function get_avg_power(solution_result::Tuple{DataFrame, Dict{Symbol, Any}})
  (df, opts) = solution_result
  if size(df, 1) == 0
    NaN
  else
    c = opts[:config]
    mean(df.power) * c.n
  end
end

function signal_sum_of_kites(values, n::Integer)
  l = length(values)
  slice = div(l, n)
  [sum([values[(i + slice * j - 1) % l + 1] for j=1:n]) for i=1:l]
end


# An angle in radians that is zero when the bridle force is pointing directly
# towards the center of rotation, and pi/2 when pointing towards the direction
# of travel
function get_bridle_alpha_b_and_force(solution_result::Tuple{DataFrame, Dict{Symbol, Any}})
end



function max_alpha_b(n_kites::Integer)
  (n_kites - 1) / n_kites * pi / 2
end
