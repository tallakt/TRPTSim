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
