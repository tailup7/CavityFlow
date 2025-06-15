import numpy as np

LX    = 1       # region length
LY    = 1
NX    = 64      # num of grids
NY    = 64
DT    = 0.00001 # time step
T_END = 1.0
U     = 10       # velocity of moving wall
NU    = 1       # kinematic viscosity
DX    = LX / (NX-1)
DY    = LY / (NY-1)
NUM_OF_PARTICLE = 1000

# coefficient matrix of Poisson equation ( A x = b )
A     = np.zeros((NX*NY, NX*NY)) 
for i in range (NX):
    for j in range (NY):
        PP = i + j * NX
        SS = i + (j - 1) * NX
        WW = (i - 1) + j * NX
        EE = (i + 1) + j * NX
        NN = i + (j + 1) * NX
        if j >=1:
            A[PP][SS] = 1 / (DY * DY)
            A[PP][PP] -= 1 / (DY * DY)
        if i >= 1:
            A[PP][WW] = 1 / (DX * DX)
            A[PP][PP] -= 1 / (DX * DX)
        if i <= NX - 2:
            A[PP][EE] = 1 / (DX * DX)
            A[PP][PP] -= 1 / (DX * DX)
        if j <= NY - 2:
            A[PP][NN] = 1 / (DY * DY)
            A[PP][PP] -= 1 / (DY * DY)