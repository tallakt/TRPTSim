
DENSITY_OF_UHWMPE_TETHER = 711.5

function tether_mass(diameter, length = 1.0)
	(diameter / 2)^2 * pi * DENSITY_OF_UHWMPE_TETHER * length
end

function tether_strength(diameter)
	diameter ^ 2 / 0.004^2 * 21_000 # using DSK90 4 mm as reference 
end

