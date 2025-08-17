import numpy as np
import config

NX    = config.NX    # num of grids
NY    = config.NY
DX    = config.DX
DY    = config.DY
DT    = config.DT

def grad(data):
    for i in range (NX):
        for j in range(NY):
            # pressure gradient term in NSeq of x direction
            # dp/dx
            if i <= NX-2:
                data.dpdx[i][j] = (data.p[i+1][j] - data.p[i][j]) / DX
            else:
                data.dpdx[i][j] = (          - data.p[i][j]) / DX
            # pressure gradient term in NSeq of y direction
            # dp/dy
            if j != NY-1:
                data.dpdy[i][j] = (data.p[i][j+1] - data.p[i][j]) / DY
            else:
                data.dpdy[i][j] = (          - data.p[i][j]) / DY

def solve_poisson(up, vp):
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
    # Right hand side of Poisson equation : b
    b    = np.zeros(NX*NY)
    for i in range (1,NX-1):
        for j in range (1,NY-1):
            b[i+j*NX] = ((-up[i-1][j]+up[i][j])/DX + (-vp[i][j-1]+vp[i][j])/DY) / DT  
    for i in range (1,NX-1):
        b[i]           = ((-up[i-1][0]    + up[i][0])/DX    + vp[i][0]   /DY) / DT
        b[i+(NY-1)*NX] = ((-up[i-1][NY-1] + up[i][NY-1])/DX - vp[i][NY-2]/DY) / DT 
    for j in range (1,NY-1):
        b[j*NX]        = (up[0][j]/DX        + (-vp[0][j-1]    + vp[0][j])   /DY) / DT
        b[NX-1+j*NX]   = (-up[NX-2][j]/DX    + (-vp[NX-1][j-1] + vp[NX-1][j])/DY) / DT 
    b[0]               = (up[0][0] / DX      + vp[0][0]      /DY) / DT   
    b[(NY-1)*NX]       = (up[0][NY-1]/DX     - vp[0][NY-2]   /DY) / DT       
    b[NX-1]            = (-up[NX-2][0]/DX    + vp[NX-1][0]   /DY) / DT          
    b[NX-1+(NY-1)*NX]  = (-up[NX-2][NY-1]/DX - vp[NX-1][NY-2]/DY) / DT
    # solve poisson eq
    phi = np.linalg.solve(A,b)
    return phi

