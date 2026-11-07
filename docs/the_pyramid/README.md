# The Pyramid: calculations

The calculations of `../the_pyramid.pdf` (AWEC 2021) were first made in a worksheet. They have
moved to samara, a computer algebra system written for this kind of work, and the samara files in
this folder are now the calculations of the paper.

| file | what it is |
|---|---|
| `pyramid.samara` | the calculation: the shaft and its moment (the MTR value), the force balance of the kite and the harvesting force for the most power, the roll controller, and the drag of the tethers with the longest shaft. It runs in about a second |
| `checks.samara` | the results of `pyramid.samara` checked against the paper and against the values of the earlier version, each line a `test_equal` that should come out true |
| `figures/*.jl` | the Julia scripts that draw the paper's plots, unchanged |

Run them with `isamara checks.samara` (it includes `pyramid.samara`).

## What the move showed

The errors found in the 2022 text and listings are in `../the_pyramid_addendum.pdf`. In short: two
values of the table of shaft geometries are truncated (C and D); the power plots and the
numerical column of the glide-number table are made at 60 degrees where the text says 30; the
"maximum" line of the shaft-length plots is the value at a twist of 90 degrees, not the maximum;
and a call in the two-section Julia listing misses an argument, hidden by a `try`/`catch`.

Two things are new beside the paper: with the elevation angle kept symbolic, the optimal MTR
value is `1 / (2 G)` for every elevation, and the optimal power is Loyd's crosswind power with
`w cos(theta_e)` for the wind; and the exact largest MTR value over the twist is
`(sqrt(l^2 - (r1 - r0)^2) - sqrt(l^2 - (r1 + r0)^2)) / (2 r1)`.

Not checked: the tether drag formulas and the shaft-length results, which rest on the companion
note on the drag of a tether with a belly.

## Figures

The seven figures of the paper are drawn by the Julia scripts in `figures/`, and a few by hand in
Sketch. samara has no plots yet, so the figures are not redrawn here: `pyramid.samara` makes
tables of the same numbers. The scripts have three corrections since 2022, so that a figure drawn
now differs from the one in the paper:

- `optimal_force.jl` draws its plots at the elevation of 30 degrees the text says, set in
  `DEFAULT_ELEVATION_ANGLE` (the 2022 figures used 60 degrees: set it to `pi / 3` to get them back).
- `zero_drag_kite_shaft_length.jl` draws the exact largest moment-to-tension-to-radius ratio
  (`exact_lambda_max`) where it drew the value at a twist of 90 degrees.
- `double_section.jl` no longer searches for buckling. In the 2022 code the search never ran (a
  call missed an argument and a `try` hid it). With the call fixed it ran, as an unbracketed secant
  search over all twists, and found twists of 180 to 3060 degrees for l0 between about 11 and 16 m
  (15/5 among them), where the plot then dropped to 0 or worse; the real buckling, below l0 = 8 m,
  it never found. The script now works the shaft out at each twist from 0 to 180 degrees: two
  sections while the tethers pull the bridle joint outwards (with the horizontal part of their
  tension, where the old check used the whole tension), one straight section r0 -> r2 once the
  bridle goes slack. The curve is within 0.1 % of the 2022 figure for l0 of 3.5 m and more, and up
  to 3.3 % higher below, where the 2022 code left out the buckled, single-section shaft.

The addendum describes why. The figures in the paper are the old ones.

## The PDF

`../the_pyramid.pdf` is the 2022 version and is not changed by this commit.
