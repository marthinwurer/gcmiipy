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

2. Set initial conditions and boundary conditions

3. do timestep and return new state

4. log state if needed

5. go to 3 if not in final state


So, we need initial conditions, a final state test, a state definition, boundary conditions, and logging.

State definition is straightforward

State is anything that changes over a timestep that cannot be easily recomputed from the previous state.

We need global constants and model constants. Global constants are based on physical laws. Model constants are constant for the model (like grid, dimensions, geography, forcings, etc).

Our initial problem is discretizing the Navier-Stokes equations.

https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations#Compressible_flow

Parameters of this formula are: 
* velocity
* density
* pressure
* change in time (in the Dt term)
* change in dimensions (for gradient terms)
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








