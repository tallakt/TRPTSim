# The Pyramid: calculations

The calculations of `../the_pyramid.pdf` (AWEC 2021) were first written as a Maple worksheet. They
have moved to samara, a computer algebra system written for
this kind of work, and the worksheet is no longer part of the repository. It is still in the git
history, as the earlier commits show.

| file | what it is |
|---|---|
| `pyramid.samara` | the calculation: the shaft and its moment (the MTR value), the force balance of the kite and the harvesting force for the most power, the roll controller, and the drag of the tethers with the longest shaft. It runs in about a second |
| `checks.samara` | the results of `pyramid.samara` checked against the paper and against the values of the earlier version, each line a `test_equal` that should come out true |
| `figures/*.jl` | the Julia scripts that draw the paper's plots, unchanged |

Run them with `isamara checks.samara` (it includes `pyramid.samara`).

## What the move showed

The samara version was written from the paper and the earlier worksheet, and checked against both.
These differences were found. Items 1 to 5 are reproduced by `pyramid.samara`; the last two are in
the Julia plot code and have not been looked at there yet.

1. **Worksheet output.** The saved output of the earlier worksheet had `tau := 0`, from a run made
   before the linear algebra package was loaded. Equations (8) and (9) of the paper are right.
2. **Generalised elevation angle.** With the elevation angle `theta_e` kept symbolic, the optimal
   power is exactly the crosswind power of Loyd with `w cos(theta_e)` for the wind, and the optimal
   MTR value is `1 / (2 G)` for every `theta_e`. This is new beside the paper.
3. **The most the shaft can carry.** Over the twist, the largest MTR value is
   `(sqrt(l^2 - (r1 - r0)^2) - sqrt(l^2 - (r1 + r0)^2)) / (2 r1)`, which goes to `r0 / l` for a
   long shaft.
4. **Table.** The values for shafts C and D are truncated, not rounded: 0.104 (the paper has 0.103)
   and 0.0578 (0.057).
5. **Tether drag.** The factor `(r0^3 + r1^3) / ((r1 - r0) r1^2)` has a sign slip against the line
   before it in the paper, and both lines leave out the moment arm of the inner cone. The factor
   is `(r0 + r1)(r0^2 + r1^2) / r1^3`, which is 4 for parallel tethers where the paper's grows
   without bound. At K = 5 the two are close (1.248 and 1.26).
6. **Shaft length (Julia).** The paper's formula for l/d has `6 pi mu`; the plot code seems to use
   `4 pi`, which gives lengths two thirds of the formula's. In the text, `mu = 0.5 MPa` should be
   GPa.
7. **Numeric MTR column (Julia).** The paper's numeric column seems read off plots cut at an MTR
   value of 0.04, and the plot code runs at an elevation of 60 degrees where the text says 30.
   Computed exactly, equation (27) holds better than the column suggests.

## Figures

The seven figures of the paper are drawn by the Julia scripts in `figures/`, and a few by hand in
Sketch. samara has no plots yet, so the figures are not redrawn here: `pyramid.samara` makes
tables of the same numbers. The scripts are left as they were, so the figures can be drawn again
as before. Where they differ from the numbers above (items 6 and 7), the figures in the paper are
the old ones.

## The PDF

`../the_pyramid.pdf` is the 2022 version and is not changed by this commit.
