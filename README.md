# GCM II Py

A recreation of GISS's GCMII Global Climate model in python

## Notes

The initial problem is to solve the Navier-Stokes equations

https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations#Compressible_flow

so, one way to do it would be to write a solver

Given an Abstract Syntax Tree (AST) of the given equations, it would solve for the
dependent variable and then automatically discretize the equation.

What would need to be done:
* Implement a solver that would do the algebra to move things around to solve 
    for the dependent variable
* implement each operation as a class
* implement each data type as a class that worked off of numpy and when 
    operations are done on them build a new AST from that

Solver would need to define discretization methods for each variable.

Like, for simplified t1 with just pgf and material derivative, you need to discretize the fields, and you need to discretize the time. 

you want to use a different method for the time and the fields. 
see http://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/01_Step_1.ipynb
actually see http://nbviewer.jupyter.org/github/barbagroup/CFDPython/blob/master/lessons/14_Step_11.ipynb
(all of them are from https://github.com/barbagroup/CFDPython)


Or I could just not do that

Manually discretize each function

Then do all the timestepping stuff too.

This would get faster results, but future stuff would be harder.

I was trying the first, but I think I want to try the second so I can get some
results.

General architecture for any of these problems:

1. define formulas for timestep

2. Define needed states, required global constants, and model constants.

3. find all PDEs and determine how they will be discretized.

4. Set initial conditions and boundary conditions

5. do timestep and return new state

6. log state if needed

7. go to 3 if not in final state


So, we need initial conditions, a final state test, a state definition, boundary conditions, and logging.

State definition is straightforward

State is anything that changes over a timestep that cannot be easily recomputed from the previous state.
(basically it's a partial with respect to time)

We need global constants and model constants. Global constants are based on physical laws. Model constants are constant for the model (like grid, dimensions, geography, forcings, etc).

Our initial problem is discretizing the Navier-Stokes equations.

https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations#Compressible_flow

Parameters of this formula are: 
* velocity
* density
* pressure
* change in time (in the Dt term)
* change in space (for gradient terms)
* viscosity
* gravity

Viscosity and gravity are global constants. They do not change between runs of the simulation.
For our model we will set viscosity and gravity to zero for the initial attempt at the 1d soultion.

The infintesimals dt and dx will be model constants, based on the timestepping and griding.

This leaves us with density, pressure, and velocity as our state.

https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)#Euler_equations

Actually we're basically solving the euler equations

We want to start with moving a tracer so we don't have to deal with changing pressure or velocity. 
We can add water vapor as a tracer, to be advected with the velocity.
That means we also need the advection equations.
https://en.wikipedia.org/wiki/Advection#The_advection_equation

Advection equations are part of the material derivative.

Finite differences will be used to solve the equations at first

https://en.wikipedia.org/wiki/Finite_difference

change in time will be forward difference

for first try, change in space will be forward difference too.

forward difference: change every dA/dx into (A(x+h) - A(x)) / (h) where h is the change in x and A is the function

backwards difference: change every dA/dx into (A(x) - A(x-h)) / (h) where h is the change in x and A is the function

central difference: change every dA/dx into (A(x+h) - A(x-h)) / (2 * h) where h is the change in x and A is the function


Do that with all partial derivatives, including gradient.




https://www.myroms.org/wiki/Time-stepping_Schemes_Review





# Stability

https://www.fluvial.ch/m/Hydraulics2_ShallowWater_2008.pdf

I'm using an explict scheme, so I have to respect the CFL conditions.

That pdf shows that the CFL for shallow water equations is:

```
CFL = (|u| + c) * dt / dx
```

Where `u` is flow velocity and `c` is wave velocity (celerity).

Shallow water is where L/h of waves is > 20.
L = wave length and h = water depth.
For shallow waves, 

```
c = sqrt(g * h)
```

For my initial conditions, h ended up being around 8029m (8km).
My minimum wavelength from dx was 300km.
```
300 / 8 = 37.5
```
I was fine being in the shallow regime.

My celerity was:
```
sqrt(G * 8km) = 280.0 m/s
```

Plug that into the CFL Formula:
```
CFL = (280 m/s) * 900s / 300km = 0.84
```

That means that my maximum `u` could be:

```
1 = (u + 280 m/s) * 900s / 300km

300km = 900u s + 252000 m

48 km = 900u s

53.3 m/s = u
```

Hmm, my `u` was set to 1 m/s. I wonder what happens if I half it?

Not enough.

https://www.io-warnemuende.de/tl_files/staff/burchard/pdf/Numerics_Shallow_Water.pdf

I found this set of notes. They told me some super important things, like all
the sets of CFL conditions for the 2d shallow water equations for all sorts of 
grids for both static and rotating flows.

The key piece of knowledge here is that the CFL needs to be less than `sqrt(1/2)`, not 1.
That comes out to be something like 0.707.
Now, if I solve for the max timestep with 1 m/s winds:

```
CFL = sqrt(1/2) <= (|u| + sqrt(g * h)) * dt / dx

sqrt(1/2) <= (1 m/s + sqrt(g * 8km)) * dt / 300km

0.707 <= (1 m/s + 280 m/s) * dt / 300km

212.1 km <= (281 m/s) * dt 

754.8 s = dt
```

I'll try 700 s.

It worked!

## momentum terms

The momentum terms actually need you to use the momentum.
This means scaling everything by the pressure/density and area.
Make sure to use the conservation form of the equations with p in the derivative.






