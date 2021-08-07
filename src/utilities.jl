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
