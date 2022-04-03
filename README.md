# TRPT Simulator

![The Pyramid Illustration](images/illustration.jpg "The Pyramid Illustration")

This is a simulator for a particular kind of TRPT type. TRPT airborne wind
energy is a class of windmills that are airborne and where the power of the
flying part of the windmill is transferred using torque transfer through a
tensile shaft. It means Tensile Rotary Power Transmission and was [described
first by Oliver
Tulloch](https://pureportal.strath.ac.uk/en/publications/tensile-rotary-power-transmission-model-development-for-airborne-)

This particular class of TRPT consists of a soft shaft (possibly with
compression rods) and a number of kites rotating in a single layer, connected
through a bridle, flying in a ring. And all the kites are of the same type.

To say what it isn't, no lifter kite is assumed, and the shaft is assumed to be
mostly pure tether and not very long. Also there is only one layer, so not
lower layers to expand the rotary shaft etc.

Why? Because these rigs are a simplest structure possible if the kite is very
controllable. So some kind of computer control of the kite during the loop is
assumed.

The simulator is intended for use with a stationary power plant as well as a
power source on a moving vehicle.


## Description of the simulation performed

The simulator is intended for quick prototyping of TRPT rigs without doing a
complete detailed simulation. So

- The shape of the layer is assumes constant and all bridle lines tensioned. It
  is a goal of the simulator though, to verify that all bridle lines are
  tensioned
- Tether drag is using a simplified model disregarding dynamics
- The roll of the kites is set directly, dynamics in roll are not accounted for
- Pitch is not modelled, it is assumed that we can set pitch to achieve the
  desired lift coeficcient as needed


## Install

In Julia, install from git with

```julia
pkg> add /path/to/TRPTSim#main
```

Then use with 

```julia
using TRPTSim
```

## Examples

First create a model/configuration

```julia
> cfg1 = ampyx_ap2()
> cfg2 = delft_lei_v3()
> cfg3 = eijkelhof()
> cfg4 = october_kite()
> cfg5 = wi_rigid_daisy()
```


Don't expect any detailed simulation. The models are coarse approximation based
on available data. And the kites are not optimized for this kind of TRPT
simulator.


You can convert a model to a `Dict` to allow Julia to pretty print the contents

```julia
> cfg |> config_to_dict
Dict{Symbol, Any} with 13 entries:
  :n              => 3
  :rho            => 1.225
  :gravity        => 9.81
  :elev           => 0.523599
  :l              => 6.7
  :safety_factor  => 3.0
  :c_d_tether     => 1.1
  :d              => 0.0012
  :s              => 0.2
  :design_c_l     => 0.635
  :c_d_fun_coeffs => [0.0148208, -0.0188157, 0.0381036, 0.0754882, 0.145568, -0.261428, -0.070654, 0.20734, -0.0597481]
  :m              => 0.42
  :radius         => 4.29
```

These values are what defines the kite parameters. To create a config with some
different parameters, do

```julia
> cfg = config(wi_rigid_daisy(), n = 6, l = 15.0)
```

Or modify a confuration by scaling it

```julia
> cfg = scale_to_area(wi_rigid_daisy(), 40.0)
```

The basic unit of simulation is the `solve` function. There is a second version
that is mostly more useful that wraps the result of the `solve` operation into
a `DataSet` for easier analysis.


```julia
> wind_speed = 12.0
> psi = 0.0
> mtr = 0.1 # moment per tension per looping radius
> force_h = 0.0
> solution = solve_sector(cfg, wind_speed, psi, mtr, force_h)
> (dataframe, solver_input) = solve_sector_df(cfg, wind_speed, psi, mtr, force_h)
```

The inputs to the solver is some parameters that are not spesifically
configuration, but more things that would notmally change. Wind speed is a good
example. The inputs besides the configuration are

- `wind_speed`: The wind speed in m/s
- `psi`: the offset in azimuth relative to pointing the AWE rig directly downwind
- `mtr`: The moment that is be transferred for a given tension
  of the shaft, divided by the looping radius. Detailed explanation below.
- `force_h`: the force generated by the implemented algorithm in a direction
  horizontal in the direction perpendicular to the shaft centerline.


The parameter `mtr` requires more explanation. It represents the shafts ability
to transfer torque given a certain tension. It is also divided by the looping
radius to make the value constant during scaling, so that a `mtr` factor is
similar across a wide range of designs and sizes.

To arrive at the actual moment transmitted through the shaft, use the formula:

> `moment = MTR * looping_radius * shaft_tension`

The `mtr` factor also represents the induction factor of the rig. An `mtr`
value of zero represents a shaft that is not transmitting any torque and the
blades are rotating unloaded. Too high of an `mtr` factor and the rig will not
perform well as the blades are slowed down.


Deciding the optimum `mtr` and tether tension is a hard optimization problem
that takes a while to calculate. Luckily, the `mtr` value will not change much,
so in this software, most times the `mtr` is just expected to be found by trial
and error, and only the tension is optimized automatically every time.

The maximum possible `mtr` value depends on the shaft geometry. It may be
calculated by calling

```julia
> radius1 = 1.0
> radius2 = 5.0
> length = 10.0
> mtr = shaft_section_mtr(radius1, radius2, length) 
```

The value of `mtr` is typically around 0.05. The maximum moment of a shoft
shaft is usually applied when the shaft section has a twist of 90 degrees. For
extending the upper wind range of a rig, the `mtr` may be increased as the
kites will need to depower to account for maximum tether loading. The other
option is to increase the tether diameter and strength. The first option
requires making the shaft more stubby or adding compressive elements, while the
latter will hurt low wind performance.


Values reported by `solve...` are usually reported per kite. For example the tension trace will report only the tension of a single kite during a full looping cycle. To get the sum of all tensions for all kites combined, use

```julia
> moment_of_all_kites = signal_sum_of_kites(dataframe.moment, solver_input[:config].n)
```

The flying speed may be estimated by this heuristic.  This is used to
initialize the solver with an initial flying speed.

```julia
> speed0 = heuristic_flying_speed(cfg, wind_speed, psi)
```

The optimal tension for a rig is estimated by iteration using. Optimal tension
is the tension that gives the maximum power output.

```julia
> tension = optimal_tension(cnf, wind_speed, psi, mtr, force_h)
```

Once you have solved the kite motion, the solution can be plotted using

```julia
> plot_solution(dataframe, solver_input)
```

A complete power curve is made and plotted like this

```julia
> winds = 8:15
> pc = power_curve(cfg, winds, psi, mtr, force_h)
> plot_power_curves(pc)
```

To plot many curves in one chart 

```julia
> plot_power_curves("First" => pc1, "Second" => pc2)
> plot_power_curves([("$mtr", power_curve(cfg, winds, psi, mtr, force_h)) for mtr = 0.01:0.01:0.08]...)
```


There is also a tension curve for mounting the windmill on a moving vessel

```julia
> plot_tension_curves(pc)
```











