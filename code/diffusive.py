import numpy as np
import config

NX    = config.NX    # num of grids
NY    = config.NY
DX    = config.DX
DY    = config.DY
DT    = config.DT
NU    = config.NU

class DiffusiveTerm:
    def __init__(self, u, v):
        self.u = u
        self.v = v
        self.diffusive_x = np.zeros([NX, NY])

    def spatial_discre_2ndorder_central(self):
        dudxe  = np.zeros([NX, NY])
        dudxw  = np.zeros([NX, NY])
        dudxdx = np.zeros([NX, NY])
        dudyn  = np.zeros([NX, NY])
        dudys  = np.zeros([NX, NY])
        dudydy = np.zeros([NX, NY])
        dvdxdx = np.zeros([NX, NY])
        dvdydy = np.zeros([NX, NY])
        dvdxe  = np.zeros([NX, NY])
        dvdxw  = np.zeros([NX, NY])
        dvdyn  = np.zeros([NX, NY])
        dvdys  = np.zeros([NX, NY])
        # diffusive term in NSeq of x direction
        for i in range (NX):
            for j in range (NY) :
                # du/dxdx
                if 1 <= i <= NX-2:
                    dudxe[i][j] = ( self.u[i+1][j] -  self.u[i][j] )   / DX 
                    dudxw[i][j] = ( self.u[i][j]   -  self.u[i-1][j] ) / DX
                elif i == 0 :
                    dudxe[i][j] = ( self.u[i+1][j] -  self.u[i][j] )   / DX
                    dudxw[i][j] = 0
                elif i == NX-1 :
                    dudxe[i][j] = 0 
                    dudxw[i][j] = ( self.u[i][j]   -  self.u[i-1][j] ) / DX
                dudxdx[i][j] = ( dudxe[i][j] - dudxw[i][j] ) / DX
                # du/dydy 
                if 1 <= j <= NY-2:
                    dudyn[i][j] = (self.u[i][j+1] - self.u[i][j]  ) / DY
                    dudys[i][j] = (self.u[i][j]   - self.u[i][j-1]) / DY
                elif j == 0:
                    dudyn[i][j] = (self.u[i][j+1] - self.u[i][j]  ) / DY
                    dudys[i][j] = 0 
                elif j == NY-1:
                    dudyn[i][j] = 0
                    dudys[i][j] = (self.u[i][j]   - self.u[i][j-1]) / DY
                dudydy[i][j] = ( dudyn[i][j] - dudys[i][j] ) / DY
        self.diffusive_x = dudxdx + dudydy
        # diffusive term in NSeq of y direction
        for i in range (NX):
            for j in range (NY) :
                # dv/dxdx
                if 1 <= i <= NX - 2:
                    dvdxe[i][j] = (self.v[i+1][j] - self.v[i][j])   / DX
                    dvdxw[i][j] = (self.v[i][j]   - self.v[i-1][j]) / DX
                elif i == 0:
                    dvdxe[i][j] = (self.v[i+1][j] - self.v[i][j])   / DX
                    dvdxw[i][j] = 0
                elif i == NX - 1:
                    dvdxe[i][j] = 0
                    dvdxw[i][j] = (self.v[i][j]   - self.v[i-1][j]) / DX
                dvdxdx[i][j] = ( dvdxe[i][j] - dvdxw[i][j] ) / DX
                # dv/dydy
                if 1 <= j <= NY - 2:
                    dvdyn[i][j] = (self.v[i][j+1] - self.v[i][j])   / DY
                    dvdys[i][j] = (self.v[i][j]   - self.v[i][j-1]) / DY
                elif j == 0 :
                    dvdyn[i][j] = (self.v[i][j+1] - self.v[i][j])   / DY
                    dvdys[i][j] = 0
                elif j == NY - 1:
                    dvdyn[i][j] = 0
                    dvdys[i][j] = (self.v[i][j]   - self.v[i][j-1]) / DY
                dvdydy[i][j] = ( dvdyn[i][j] - dvdys[i][j] ) / DY
        self.diffusive_y = dvdxdx + dvdydy

    def temporal_discre_2ndorder_adamsbashforth(self, diffusive_pre):
        diffusive_x = 1.5 * self.diffusive_x - 0.5*diffusive_pre.diffusive_x
        diffusive_y = 1.5 * self.diffusive_y - 0.5*diffusive_pre.diffusive_y
        return diffusive_x, diffusive_y
        
    def temporal_discre_3rdorder_adamsbashforth(self, diffusive_pre, diffusive_prepre):
        diffusive_x = 23/12 * self.diffusive_x - 16/12 * diffusive_pre.diffusive_x + 5/12 * diffusive_prepre.diffusive_x
        diffusive_y = 23/12 * self.diffusive_y - 16/12 * diffusive_pre.diffusive_y + 5/12 * diffusive_prepre.diffusive_y
        return diffusive_x, diffusive_y


    # def temporal_discre_2ndorder_cranknicolson():


# solve elliptic partial differential equation of x direction 
def prepare_coef_matrix_of_epde():
    coef_matrix = np.zeros((NX*NY, NX*NY)) 
    for i in range(NX):
        for j in range(NY):
            if 1 <= i <= NX-2 and 1 <= j <= NY-2:
                coef_matrix[i * NY + j][(i-1)*NY + j  ] = -0.5*NU*DT / (DX)**2
                coef_matrix[i * NY + j][ i   *NY + j-1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + NU*DT *(1/(DX)**2 + 1/(DY)**2)
                coef_matrix[i * NY + j][ i   *NY + j+1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][(i+1)*NY + j]   = -0.5*NU*DT / (DX)**2
            elif 1 <= i <= NX-2 and j == 0 :
                coef_matrix[i * NY + j][(i-1)*NY + j  ] = -0.5*NU*DT / (DX)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(2/(DX)**2 + 1/(DY)**2)
                coef_matrix[i * NY + j][ i   *NY + j+1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][(i+1)*NY + j]   = -0.5*NU*DT / (DX)**2
            elif 1 <= i <= NX-2 and j == NY-1 :
                coef_matrix[i * NY + j][(i-1)*NY + j  ] = -0.5*NU*DT / (DX)**2
                coef_matrix[i * NY + j][ i   *NY + j-1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(2/(DX)**2 + 1/(DY)**2)
                coef_matrix[i * NY + j][(i+1)*NY + j]   = -0.5*NU*DT / (DX)**2
            elif i == 0 and 1 <= j <= NY-2:
                coef_matrix[i * NY + j][ i   *NY + j-1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(1/(DX)**2 + 2/(DY)**2)
                coef_matrix[i * NY + j][ i   *NY + j+1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][(i+1)*NY + j]   = -0.5*NU*DT / (DX)**2
            elif i == NX-1 and 1 <= j <= NY-2:
                coef_matrix[i * NY + j][(i-1)*NY + j  ] = -0.5*NU*DT / (DX)**2
                coef_matrix[i * NY + j][ i   *NY + j-1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(1/(DX)**2 + 2/(DY)**2)
                coef_matrix[i * NY + j][ i   *NY + j+1] = -0.5*NU*DT / (DY)**2
            elif i == 0 and j == 0:
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(1/(DX)**2 + 1/(DY)**2)
                coef_matrix[i * NY + j][ i   *NY + j+1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][(i+1)*NY + j]   = -0.5*NU*DT / (DX)**2
            elif i == 0 and j == NY-1: 
                coef_matrix[i * NY + j][ i   *NY + j-1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(1/(DX)**2 + 1/(DY)**2)
                coef_matrix[i * NY + j][(i+1)*NY + j]   = -0.5*NU*DT / (DX)**2
            elif i == NX-1 and j == 0:
                coef_matrix[i * NY + j][(i-1)*NY + j  ] = -0.5*NU*DT / (DX)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(1/(DX)**2 + 1/(DY)**2)
                coef_matrix[i * NY + j][ i   *NY + j+1] = -0.5*NU*DT / (DY)**2
            elif i == NX-1 and j == NY-1:
                coef_matrix[i * NY + j][(i-1)*NY + j  ] = -0.5*NU*DT / (DX)**2
                coef_matrix[i * NY + j][ i   *NY + j-1] = -0.5*NU*DT / (DY)**2
                coef_matrix[i * NY + j][ i   *NY + j]   = 1 + 0.5*NU*DT *(1/(DX)**2 + 1/(DY)**2)
    return coef_matrix

