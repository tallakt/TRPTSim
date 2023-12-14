
function shaft_section_c_d(l::Number, d::Number, radius_1::Number, radius_2::Number, kite_s::Number; tether_c_d::Number = 1.1)
    tether_c_d * d * l / (3 * kite_s) * (radius_1^2 + radius_1 * radius_2 + radius_2^2) / radius_2^2
end


function shaft_section_m_t_factor(radius_1::Number, radius_2::Number, l::Number)
    radius_1 * radius_2 / sqrt(2 * radius_1 * radius_2 - radius_1^2 - radius_2^2 + l^2)
end


function shaft_section_mtr(radius_1::Number, radius_2::Number, l::Number)
    radius_1 / sqrt(l^2 - radius_1^2 - radius_2^2)
end


function shaft_section_cone_angle(radius_1::Number, radius_2::Number, l::Number)
    asin((radius_2 - radius_1) / l)
end


function shaft_section_max_length(radius_1::Number, radius_2::Number, l::Number)
    l * cos(shaft_section_cone_angle(radius_1, radius_2, l))
end


function shaft_section_compression(radius_1::Number, radius_2::Number, l::Number, delta::Number)
    # the tension of the shaft tether is shaft_tension / n_kites / compression
    # The compression is the length of the shaft vs the length of the tether for
    # a single shaft section
    l / sqrt(2 * radius_1 * radius_2 - radius_1^2 - radius_2^2 + l^2)
end
