################################################
# 2D cavity flow with finite difference method
#
# coupling : SMAC method
#
# space    : 2nd order central
#
# time     : 
#
# grid     : staggaered
#
################################################

import config
import numpy as np
import convective
import diffusive
from dataclasses import dataclass, field
import pressure

NX    = config.NX    # num of grids
NY    = config.NY
DT    = config.DT
U     = config.U     # velocity of moving wall
NU    = config.NU    # kinematic viscosity
DX    = config.DX
DY    = config.DY

@dataclass
class StepData:
    t       : float      = 0.0
    u       : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))
    v       : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))
    p       : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))
    conv    : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))
    diff    : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))
    dpdx    : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))
    dpdy    : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))
    u_magni : np.ndarray = field(default_factory=lambda: np.zeros((NX, NY)))

def set_boundary_condition(data):
    for i in range(NX):
        if 1 <= i <= NX-2:
            data.u[i][NY-1] = U
        data.u[i][0]    = 0
        data.v[i][NY-1] = 0
        data.v[i][0]    = 0
    for j in range (NY):
        data.u[0][j]    = 0
        data.v[0][j]    = 0
        data.u[NX-1][j] = 0
        data.v[NX-1][j] = 0
    return data

def calc_onestep(data_prev, data_prev2, data_prev3):
    data_current = StepData(t = data_prev.t + DT)
    # set temporary variables
    u    = np.zeros((NX, NY)) 
    v    = np.zeros((NX, NY))
    up   = np.zeros((NX, NY))    # u prediction
    vp   = np.zeros((NX, NY))    

    data_current = set_boundary_condition(data_current)

    # 
    if data_current.t == DT:
        data_prev.conv = convective.ConvectiveTerm(data_prev.u,data_prev.v)
        data_prev.conv.spatial_discre_2ndorder_central()
        convterm_x = data_prev.conv.convective_x
        convterm_y = data_prev.conv.convective_y
        data_prev.diff = diffusive.DiffusiveTerm(data_prev.u,data_prev.v)
        data_prev.diff.spatial_discre_2ndorder_central() 
        diffterm_x = data_prev.diff.diffusive_x
        diffterm_y = data_prev.diff.diffusive_y
    
    elif data_current.t == 2 * DT:
        data_prev.conv = convective.ConvectiveTerm(data_prev.u,data_prev.v)
        data_prev.conv.spatial_discre_2ndorder_central()
        convterm_x, convterm_y = data_prev.conv.temporal_discre_2ndorder_adamsbashforth(data_prev2.conv)
        data_prev.diff = diffusive.DiffusiveTerm(u,v)
        data_prev.diff.spatial_discre_2ndorder_central()
        diffterm_x, diffterm_y = data_prev.diff.temporal_discre_2ndorder_adamsbashforth(data_prev2.diff)
        
    else:
        data_prev.conv = convective.ConvectiveTerm(data_prev.u,data_prev.v)
        data_prev.conv.spatial_discre_2ndorder_central()
        convterm_x, convterm_y = data_prev.conv.temporal_discre_3rdorder_adamsbashforth(data_prev2.conv, data_prev3.conv)
        data_prev.diff = diffusive.DiffusiveTerm(u,v)
        data_prev.diff.spatial_discre_2ndorder_central()
        diffterm_x, diffterm_y = data_prev.diff.temporal_discre_3rdorder_adamsbashforth(data_prev2.diff, data_prev3.diff)

    pressure.grad(data_prev)

    # u prediction
    up = data_prev.u + DT * (- convterm_x - data_prev.dpdx + NU * diffterm_x)
    # v prediction
    vp = data_prev.v + DT * (- convterm_y - data_prev.dpdy + NU * diffterm_y) 

    # solve Poisson equation      
    phi = pressure.solve_poisson(up, vp)

    # modify u, v, p
    for i in range (NX):
        for j in range (NY):
            data_current.p[i][j] = data_prev.p[i][j] + phi[i+j*NX]   
    for i in range (1, NX-1):
        for j in range (1, NY-1):
            data_current.u[i][j] = up[i][j] - DT * (phi[i+1+j*NX]   - phi[i+j*NX]) / DX    
    for i in range (1, NX-1):
        for j in range(1, NY-1):    
            data_current.v[i][j] = vp[i][j] - DT * (phi[i+(j+1)*NX] - phi[i+j*NX]) / DY 


    # rename pre, pre2, and pre3 for next step
    if data_current.t == DT:
        data_prev2 = data_prev
        data_prev  = data_current
    else:
        data_prev3 = data_prev2
        data_prev2 = data_prev
        data_prev  = data_current

    return data_current, data_prev, data_prev2, data_prev3

def calc_onestep_implicit(data_prev, data_prev2, data_prev3):
    data_current = StepData(t = data_prev.t + DT)
    # set temporary variables
    u    = np.zeros((NX, NY)) 
    v    = np.zeros((NX, NY))
    up   = np.zeros((NX, NY))    # u prediction
    vp   = np.zeros((NX, NY))    

    data_current = set_boundary_condition(data_current)

    # 
    if data_current.t == DT:
        data_prev.conv = convective.ConvectiveTerm(data_prev.u,data_prev.v)
        data_prev.conv.spatial_discre_2ndorder_central()
        convterm_x = data_prev.conv.convective_x
        convterm_y = data_prev.conv.convective_y
        data_prev.diff = diffusive.DiffusiveTerm(data_prev.u,data_prev.v)
        data_prev.diff.spatial_discre_2ndorder_central() 
    
    elif data_current.t == 2 * DT:
        data_prev.conv = convective.ConvectiveTerm(data_prev.u,data_prev.v)
        data_prev.conv.spatial_discre_2ndorder_central()
        convterm_x, convterm_y = data_prev.conv.temporal_discre_2ndorder_adamsbashforth(data_prev2.conv)
        data_prev.diff = diffusive.DiffusiveTerm(u,v)
        data_prev.diff.spatial_discre_2ndorder_central()
        
    else:
        data_prev.conv = convective.ConvectiveTerm(data_prev.u,data_prev.v)
        data_prev.conv.spatial_discre_2ndorder_central()
        convterm_x, convterm_y = data_prev.conv.temporal_discre_3rdorder_adamsbashforth(data_prev2.conv, data_prev3.conv)
        data_prev.diff = diffusive.DiffusiveTerm(u,v)
        data_prev.diff.spatial_discre_2ndorder_central()

    pressure.grad(data_prev)

    righthand_of_epde_x = np.zeros(NX*NY)
    righthand_of_epde_y = np.zeros(NX*NY)
    coef_matrix = diffusive.prepare_coef_matrix_of_epde()
    for i in range(NX):
        for j in range(NY):
            righthand_of_epde_x[i*NY+j] = data_prev.u[i][j] + DT*(- convterm_x[i][j] - data_prev.dpdx[i][j] + 0.5* NU * data_prev.diff.diffusive_x[i][j])
            righthand_of_epde_y[i*NY+j] = data_prev.v[i][j] + DT*(- convterm_y[i][j] - data_prev.dpdy[i][j] + 0.5* NU * data_prev.diff.diffusive_y[i][j])
    up_column = np.linalg.solve(coef_matrix, righthand_of_epde_x)
    vp_column = np.linalg.solve(coef_matrix, righthand_of_epde_y)
    for i in range(NX):
        for j in  range(NY):
            up[i][j] = up_column[i*NY + j]
            vp[i][j] = vp_column[i*NY + j]

    # solve Poisson equation      
    phi = pressure.solve_poisson(up, vp)
    phi_matrix   = np.zeros((NX, NY))
    for i in range(NX):
        for j in range(NY):
            phi_matrix[i][j] = phi[i+j*NX]


    # modify u, v, p
    for i in range (NX):
        for j in range (NY):
            if 1 <= i <= NX-2 and 1 <= j <= NY-2:
                phi_laplacian = -0.5*NU*DT*(phi_matrix[i-1][j]/(DX)**2 + phi_matrix[i][j-1]/(DY)**2 -2* phi_matrix[i][j]*(1/(DX)**2 + 1/(DY)**2)\
                                            + phi_matrix[i][j+1]/(DY)**2 + phi_matrix[i+1][j]/(DX)**2 )
            elif i == 0 and 1 <= j <= NY-2:
                phi_laplacian = -0.5*NU*DT*( phi_matrix[i][j-1]/(DY)**2 - phi_matrix[i][j]*(1/(DX)**2 + 2/(DY)**2)\
                            + phi_matrix[i][j+1]/(DY)**2 + phi_matrix[i+1][j]/(DX)**2 )
            elif i == NX -1 and 1 <= j <= NY -2:
                phi_laplacian = -0.5*NU*DT*(phi_matrix[i-1][j]/(DX)**2 + phi_matrix[i][j-1]/(DY)**2 - phi_matrix[i][j]*(1/(DX)**2 + 2/(DY)**2)\
                                            + phi_matrix[i][j+1]/(DY)**2  )
            elif 1 <= i <= NX-2 and j ==0:
                phi_laplacian = -0.5*NU*DT*(phi_matrix[i-1][j]/(DX)**2 - phi_matrix[i][j]*(2/(DX)**2 + 1/(DY)**2)\
                                            + phi_matrix[i][j+1]/(DY)**2 + phi_matrix[i+1][j]/(DX)**2 )
            elif 1 <= i <= NX-2 and j == NY-1:
                phi_laplacian = -0.5*NU*DT*(phi_matrix[i-1][j]/(DX)**2 + phi_matrix[i][j-1]/(DY)**2 - phi_matrix[i][j]*(2/(DX)**2 + 1/(DY)**2)\
                                            + phi_matrix[i+1][j]/(DX)**2 )
            elif i == 0 and j == 0:
                phi_laplacian = -0.5*NU*DT*(- phi_matrix[i][j]*(1/(DX)**2 + 1/(DY)**2) + phi_matrix[i][j+1]/(DY)**2 + phi_matrix[i+1][j]/(DX)**2 )
            elif i == 0 and j == NY-1:
                phi_laplacian = -0.5*NU*DT*( phi_matrix[i][j-1]/(DY)**2 - phi_matrix[i][j]*(1/(DX)**2 + 1/(DY)**2) + phi_matrix[i+1][j]/(DX)**2 )
            elif i == NX-1 and j ==0:
                phi_laplacian = -0.5*NU*DT*(phi_matrix[i-1][j]/(DX)**2 - phi_matrix[i][j]*(1/(DX)**2 + 1/(DY)**2) + phi_matrix[i][j+1]/(DY)**2 )
            elif i == NX-1 and j == NY-1:
                phi_laplacian = -0.5*NU*DT*(phi_matrix[i-1][j]/(DX)**2 + phi_matrix[i][j-1]/(DY)**2 - phi_matrix[i][j]*(1/(DX)**2 + 1/(DY)**2))
            
            data_current.p[i][j] = data_prev.p[i][j] + phi[i+j*NX] + phi_laplacian
    for i in range (1, NX-1):
        for j in range (1, NY-1):
            data_current.u[i][j] = up[i][j] - DT * (phi[i+1+j*NX]   - phi[i+j*NX]) / DX    
    for i in range (1, NX-1):
        for j in range(1, NY-1):    
            data_current.v[i][j] = vp[i][j] - DT * (phi[i+(j+1)*NX] - phi[i+j*NX]) / DY 
    for i in range(NX):
        for j in range(NY):
            data_current.u_magni[i][j] = np.sqrt(data_current.u[i][j]**2 + data_current.v[i][j]**2)

    # rename pre, pre2, and pre3 for next step
    if data_current.t == DT:
        data_prev2 = data_prev
        data_prev  = data_current
    else:
        data_prev3 = data_prev2
        data_prev2 = data_prev
        data_prev  = data_current

    return data_current, data_prev, data_prev2, data_prev3