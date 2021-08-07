
function nan_to_value(x, value)
  isnan(x) ? value : x
end


function grid_optimize_1d_helper(xs, ys, fun::Function, iterations::Integer)::Float64
  if(iterations == 0)
    (_, i) = findmin(nan_to_value.(ys, Inf))
    xs[i]
  else
    xs_ = [xs[1], (xs[1] + xs[2]) * 0.5, xs[2], (xs[2] + xs[3]) * 0.5, xs[3]]
    ys_ = [ys[1], fun(xs_[2]), ys[2], fun(xs_[4]), ys[3]]
    (_, i) = findmin(nan_to_value.(ys_, Inf))
    i_ = clamp((1:length(ys_))[i], 2, length(ys_) - 1)
    grid_optimize_1d_helper(xs_[(i_-1):(i_+1)], ys_[(i_-1):(i_+1)], fun, iterations - 1)
  end
end


function grid_optimize_1d(xs, fun::Function; iterations::Integer = 10)::Float64
  ys = [fun(x) for x=xs]
  if all(isnan, ys)
    NaN
  else
    (_, i) = findmin(nan_to_value.(ys, Inf))
    i_ = clamp((1:length(xs))[i], 2, length(ys) - 1)
    grid_optimize_1d_helper(xs[(i_-1):(i_+1)], ys[(i_-1):(i_+1)], fun, iterations)
  end

end


