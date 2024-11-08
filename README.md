# Cavity Flow

## Overview
Cavity flow is a flow in which a horizontal flow crosses the upper side, for example, a hollow in the bottom of a river.
There is no need to consider inflow and outflow conditions and the region is simple in shape, 
so it is well known as a benchmark problem in fluid analysis.

## Governing Equation
![bbb](https://github.com/user-attachments/assets/827f069f-1aa8-45b1-be93-62e2700f3ffd)
<br>
boundary condition : u = v = 0 at east, west, south wall. u = U, v = 0 at south wall. <br>
initial condition : uniformly 0

## Result
![cavity](https://github.com/user-attachments/assets/0bbca827-3d93-4880-bc35-d0555a9d3053)
<br>
<br>
In this calculation, U (which is the velocity of upper side) is set 1, and Re = 1.<br>
So, the flow is very slow, and the flow field reaches a steady state.

### Comment
In this program, 1st order Euler explicit method is employed about time to make it easier to code. <br>
Therefore, the program can analyze only flows with very low Reynolds numbers. If U is set a little higher or the number of grids is increased, the calculation diverges. <br> 
So, the program will be modified to Euler implicit solution. Also, it will be added code of particle tracking to make the flow field more easily visualized .
