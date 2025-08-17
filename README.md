# Cavity Flow

## Overview
Cavity flow is a flow in which a horizontal flow crosses the upper side, for example, a hollow in the bottom of a river.
There is no need to consider inflow and outflow conditions and the region is simple in shape, 
so it is well known as a benchmark problem in fluid analysis.

## Governing Equation
![bbb](https://github.com/user-attachments/assets/827f069f-1aa8-45b1-be93-62e2700f3ffd)
<br>
boundary condition : u = v = 0 at east, west, south wall. u = U, v = 0 at north wall. <br>
initial condition : uniformly 0

## Result
![cavity](https://github.com/user-attachments/assets/0bbca827-3d93-4880-bc35-d0555a9d3053)
<br>
<br>
In this calculation, U (which is the velocity of upper side, as a reference velocity) is 50, and Re = 50.<br>
So, the flow is slow, and the flow field reaches a steady state.


### Particle Tracking
![ani](https://github.com/user-attachments/assets/5c09915e-432e-43db-9a99-b1fe968eb33b)
Particles were scattered randamly to make the flow easier to see (ex, recirculation region). <br>
These particles have no mass and do not affect the flow field.
