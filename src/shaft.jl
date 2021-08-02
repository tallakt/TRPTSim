
function shaft_c_d(l, d, radius_1, radius_2, kite_s; tether_c_d = 1.1)
  tether_c_d * d * l / (3 * kite_s) * (radius_1^2 + radius_1 * radius_2 + radius_2^2) / radius_2^2
end


function shaft_m_t_factor(radius_1, radius_2, l)
  radius_1 * radius_2 / sqrt(2 * radius_1 * radius_2 - radius_1^2 - radius_2^2 + l^2)
end
