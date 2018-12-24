# GCM II Py

A recreation of GISS's GCMII Global Climate model in python

## Notes

The end is near

so, one way to do it would be to write a solver

Given an Abstract Syntax Tree (AST) of the given equations, it would solve for the
dependent variable and then automatically discretize the equation.

What would need to be done:
* Implement a solver that would do the algebra to move things around to solve 
    for the dependent variable
* implement each operation as a class
* implement each data type as a class that worked off of numpy and when 
    operations are done on them build a new AST from that

Or I could just not do that

Manually discretize each function

Then do all the timestepping stuff too.

This would get faster results, but future stuff would be harder.

I was trying the first, but I think I want to try the second so I can get some
results.