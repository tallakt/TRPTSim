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
- `double_section.jl` calls `lambda` with its fourth argument in the search for buckling, and an
  error in the code is no longer hidden by the `try`/`catch`. For r0 = r1 = 1, r2 = 4 and the
  splits 10/10, 5/15 and 15/5 the search finds no buckling, and the results are as before.

The addendum describes why. The figures in the paper are the old ones.

## The PDF

`../the_pyramid.pdf` is the 2022 version and is not changed by this commit.
