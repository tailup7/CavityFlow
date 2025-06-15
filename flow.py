################################################
# 2D cavity flow with finite difference method
#
# coupling : SMAC method
#
# space    : 2nd order central
#
# time     : 1st order Euler   # TODO:陰解法にしてレイノルズ数の高い計算を行う
#
# grid     : staggaered
#
################################################

import config
import numpy as np

NX    = config.NX    # num of grids
NY    = config.NY
DT    = config.DT
U     = config.U     # velocity of moving wall
NU    = config.NU    # kinematic viscosity
DX    = config.DX
DY    = config.DY
A     = config.A     # coefficient matrix of Poisson equation ( A x = b )

def calc_onestep(u, v, p):
    
    # set temporary variables
    up   = np.zeros((NX, NY))    # u prediction
    vp   = np.zeros((NX, NY))  
    b    = np.zeros(NX*NY)       # righthand side of Poisson equation ( A x = b )
    dpdx = np.zeros((NX, NY))  
    dpdy = np.zeros((NX, NY))   
    phi  = np.zeros(NX*NY)       # pressure modification 

    # set boundary condition
    for i in range(NX):
        if 1 <= i <= NX-2:
            u[i][NY-1] = U
        u[i][0]    = 0
        v[i][NY-1] = 0
        v[i][0]    = 0
    for j in range (NY):
        u[0][j]    = 0
        v[0][j]    = 0
        u[NX-1][j] = 0
        v[NX-1][j] = 0

    for i in range (NX):
        for j in range (NY) :
            # convection term in NSeq of x direction
            # d(uu)/dx
            if 1 <= i <= NX-2:
                ue = (u[i+1][j] + u[i][j])   / 2
                uw = (u[i][j]   + u[i-1][j]) / 2 
            elif i == 0:
                ue = (u[i+1][j] + u[i][j])   / 2
                uw = 0 
            elif i == NX-1:
                ue = 0
                uw = (u[i][j]   + u[i-1][j]) / 2 
            duudx = (ue**2 - uw **2) / DX
            # d(uv)/dy
            if i != NX-1:
                if 1 <= j <= NY-2:
                    un = ( u[i][j]   + u[i][j+1] ) / 2
                    us = ( u[i][j]   + u[i][j-1] ) / 2
                    vn = ( v[i][j]   + v[i+1][j] ) / 2
                    vs = ( v[i][j-1] + v[i+1][j-1] ) / 2
                elif j == 0 :
                    un = ( u[i][j]   + u[i][j+1] ) / 2
                    us = 0
                    vn = ( v[i][j]   + v[i+1][j] ) / 2
                    vs = 0
                elif j == NY-1:
                    un = ( u[i][j]               ) / 2
                    us = ( u[i][j]   + u[i][j-1] ) / 2
                    vn = ( v[i][j]   + v[i+1][j] ) / 2
                    vs = ( v[i][j-1] + v[i+1][j-1] ) / 2
            else:
                un = 0
                us = 0
                vn = 0
                vs = 0 
            duvdy = (un*vn - us*vs) / DY

            # pressure gradient term in NSeq of x direction
            # dp/dx
            if i <= NX-2:
                dpdx[i][j] = (p[i+1][j] - p[i][j]) / DX
            else:
                dpdx[i][j] = (          - p[i][j]) / DX

            # diffusive term in NSeq of x direction
            # du/dxdx
            if 1 <= i <= NX-2:
                dudxe = ( u[i+1][j] -  u[i][j] )   / DX 
                dudxw = ( u[i][j]   -  u[i-1][j] ) / DX
            elif i == 0 :
                dudxe = ( u[i+1][j] -  u[i][j] )   / DX
                dudxw = 0
            elif i == NX-1 :
                dudxe = 0 
                dudxw = ( u[i][j]   -  u[i-1][j] ) / DX
            dudxdx = ( dudxe - dudxw ) / DX
            # du/dydy 
            if 1 <= j <= NY-2:
                dudyn = (u[i][j+1] - u[i][j]  ) / DY
                dudys = (u[i][j]   - u[i][j-1]) / DY
            elif j == 0:
                dudyn = (u[i][j+1] - u[i][j]  ) / DY
                dudys = 0 
            elif j == NY-1:
                dudyn = 0
                dudys = (u[i][j]   - u[i][j-1]) / DY
            dudydy = ( dudyn - dudys ) / DY

            # u prediction
            up[i][j] = u[i][j] + DT * (- (duudx + duvdy) - dpdx[i][j] + NU * (dudxdx + dudydy))

            # convection term in NSeq of y direction
            # d(uv)/dx
            if j != NY-1:
                if 1 <= i <= NX - 2:
                    ue = (u[i][j]   + u[i][j+1]  ) / 2
                    uw = (u[i-1][j] + u[i-1][j+1]) / 2
                    ve = (v[i][j]   + v[i+1][j]  ) / 2
                    vw = (v[i][j]   + v[i-1][j]  ) / 2
                elif i == 0:
                    ue = (u[i][j]   + u[i][j+1]  ) / 2
                    uw = 0
                    ve = (v[i][j]   + v[i+1][j]  ) / 2
                    vw = 0
                elif i == NX - 1:
                    ue = 0
                    uw = (u[i-1][j] + u[i-1][j+1]) / 2
                    ve = 0
                    vw = (v[i][j]   + v[i-1][j]  ) / 2
            else:
                if i != 0:
                    ue = u[i][j]
                    uw = u[i-1][j]
                    ve = 0 
                    vw = 0
                if i == 0:
                    ue = u[i][j]
                    uw = 0
                    ve = 0 
                    vw = 0
            duvdx = ( ue * ve - uw * vw)
            # d(vv)/dy
            if 1 <= j <= NY - 2:
                vn = (v[i][j] + v[i][j+1]) / 2
                vs = (v[i][j] + v[i][j-1]) / 2
            elif j == 0:
                vn = (v[i][j] + v[i][j+1]) / 2
                vs = 0
            elif j == NY - 1:
                vn = 0 
                vs = (v[i][j] + v[i][j-1]) / 2
            dvvdy = ( vn*vn - vs*vs ) / DY

            # pressure gradient term in NSeq of y direction
            # dp/dy
            if j != NY-1:
                dpdy[i][j] = (p[i][j+1] - p[i][j]) / DY
            else:
                dpdy[i][j] = (          - p[i][j]) / DY

            # diffusive term in NSeq of y direction
            # dv/dxdx
            if 1 <= i <= NX - 2:
                dvdxe = (v[i+1][j] - v[i][j])   / DX
                dvdxw = (v[i][j]   - v[i-1][j]) / DX
            elif i == 0:
                dvdxe = (v[i+1][j] - v[i][j])   / DX
                dvdxw = 0
            elif i == NX - 1:
                dvdxe = 0
                dvdxw = (v[i][j]   - v[i-1][j]) / DX
            dvdxdx = ( dvdxe - dvdxw ) / DX
            # dv/dydy
            if 1 <= j <= NY - 2:
                dvdyn = (v[i][j+1] - v[i][j])   / DY
                dvdys = (v[i][j]   - v[i][j-1]) / DY
            elif j == 0 :
                dvdyn = (v[i][j+1] - v[i][j])   / DY
                dvdys = 0
            elif j == NY - 1:
                dvdyn = 0
                dvdys = (v[i][j]   - v[i][j-1]) / DY
            dvdydy = ( dvdyn - dvdys ) / DY

            # v prediction
            vp[i][j] = v[i][j] + DT * ( - (duvdx + dvvdy) - dpdy[i][j] + NU * (dvdxdx + dvdydy)) 

    # Right side of Poisson equation : b
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

    # solve Poisson equation      
    phi = np.linalg.solve(A,b)

    # modify u, v, p
    for i in range (NX):
        for j in range (NY):
            p[i][j] += phi[i+j*NX]   
    for i in range (NX-1):
        for j in range (NY):
            u[i][j] = up[i][j] - DT * (phi[i+1+j*NX]   - phi[i+j*NX]) / DX    
    for i in range (NX):
        for j in range(NY-1):
            v[i][j] = vp[i][j] - DT * (phi[i+(j+1)*NX] - phi[i+j*NX]) / DY 

    return u, v, p