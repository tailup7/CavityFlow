import numpy as np
import config

NX    = config.NX    # num of grids
NY    = config.NY
DX    = config.DX
DY    = config.DY

class ConvectiveTerm:
    def __init__(self, u, v):
        self.u = u
        self.v = v
        self.convective_x = np.zeros([NX, NY])
        self.convective_y = np.zeros([NX, NY])

    def spatial_discre_2ndorder_central(self):
        ue    = np.zeros([NX, NY])
        uw    = np.zeros([NX, NY])
        un    = np.zeros([NX, NY])
        us    = np.zeros([NX, NY])
        vn    = np.zeros([NX, NY])
        vs    = np.zeros([NX, NY])
        ve    = np.zeros([NX, NY])
        vw    = np.zeros([NX, NY])
        duudx = np.zeros([NX, NY])
        dvvdy = np.zeros([NX, NY])
        duvdx = np.zeros([NX, NY])
        duvdy = np.zeros([NX, NY])
        # # convective term in NSeq of x direction
        for i in range (NX):
            for j in range (NY) :
                # d(uu)/dx
                if 1 <= i <= NX-2:
                    ue[i][j] = (self.u[i+1][j] + self.u[i][j])   / 2
                    uw[i][j] = (self.u[i][j]   + self.u[i-1][j]) / 2 
                elif i == 0:
                    ue[i][j] = (self.u[i+1][j] + self.u[i][j])   / 2
                    uw[i][j] = 0 
                elif i == NX-1:
                    ue[i][j] = 0
                    uw[i][j] = (self.u[i][j]   + self.u[i-1][j]) / 2 
                duudx[i][j] = (ue[i][j] **2 - uw[i][j] **2) / DX
                # d(uv)/dy
                if i != NX-1:
                    if 1 <= j <= NY-2:
                        un[i][j] = ( self.u[i][j]   + self.u[i][j+1] ) / 2
                        us[i][j] = ( self.u[i][j]   + self.u[i][j-1] ) / 2
                        vn[i][j] = ( self.v[i][j]   + self.v[i+1][j] ) / 2
                        vs[i][j] = ( self.v[i][j-1] + self.v[i+1][j-1] ) / 2
                    elif j == 0 :
                        un[i][j] = ( self.u[i][j]   + self.u[i][j+1] ) / 2
                        us[i][j] = 0
                        vn[i][j] = ( self.v[i][j]   + self.v[i+1][j] ) / 2
                        vs[i][j] = 0
                    elif j == NY-1:
                        un[i][j] = ( self.u[i][j]               ) / 2
                        us[i][j] = ( self.u[i][j]   + self.u[i][j-1] ) / 2
                        vn[i][j] = ( self.v[i][j]   + self.v[i+1][j] ) / 2
                        vs[i][j] = ( self.v[i][j-1] + self.v[i+1][j-1] ) / 2
                else:
                    un[i][j] = 0
                    us[i][j] = 0
                    vn[i][j] = 0
                    vs[i][j] = 0 
                duvdy[i][j] = (un[i][j] * vn[i][j] - us[i][j] * vs[i][j]) / DY
        self.convective_x = duudx + duvdy
        # # convective term in NSeq of y direction
        for i in range (NX):
            for j in range (NY) :
                # d(uv)/dx
                if j != NY-1:
                    if 1 <= i <= NX - 2:
                        ue[i][j] = (self.u[i][j]   + self.u[i][j+1]  ) / 2
                        uw[i][j] = (self.u[i-1][j] + self.u[i-1][j+1]) / 2
                        ve[i][j] = (self.v[i][j]   + self.v[i+1][j]  ) / 2
                        vw[i][j] = (self.v[i][j]   + self.v[i-1][j]  ) / 2
                    elif i == 0:
                        ue[i][j] = (self.u[i][j]   + self.u[i][j+1]  ) / 2
                        uw[i][j] = 0
                        ve[i][j] = (self.v[i][j]   + self.v[i+1][j]  ) / 2
                        vw[i][j] = 0
                    elif i == NX - 1:
                        ue[i][j] = 0
                        uw[i][j] = (self.u[i-1][j] + self.u[i-1][j+1]) / 2
                        ve[i][j] = 0
                        vw[i][j] = (self.v[i][j]   + self.v[i-1][j]  ) / 2
                else:
                    if i != 0:
                        ue[i][j] = self.u[i][j]
                        uw[i][j] = self.u[i-1][j]
                        ve[i][j] = 0 
                        vw[i][j] = 0
                    if i == 0:
                        ue[i][j] = self.u[i][j]
                        uw[i][j] = 0
                        ve[i][j] = 0 
                        vw[i][j] = 0
                duvdx[i][j] = ( ue[i][j] * ve[i][j] - uw[i][j] * vw[i][j])
                # d(vv)/dy
                if 1 <= j <= NY - 2:
                    vn[i][j] = (self.v[i][j] + self.v[i][j+1]) / 2
                    vs[i][j] = (self.v[i][j] + self.v[i][j-1]) / 2
                elif j == 0:
                    vn[i][j] = (self.v[i][j] + self.v[i][j+1]) / 2
                    vs[i][j] = 0
                elif j == NY - 1:
                    vn[i][j] = 0 
                    vs[i][j] = (self.v[i][j] + self.v[i][j-1]) / 2
                dvvdy[i][j] = ( vn[i][j] ** 2 - vs[i][j] ** 2 ) / DY
        self.convective_y = duvdx + dvvdy

    def temporal_discre_2ndorder_adamsbashforth(self, convective_prev):
        convective_x = 1.5 * self.convective_x - 0.5*convective_prev.convective_x
        convective_y = 1.5 * self.convective_y - 0.5*convective_prev.convective_y
        return convective_x, convective_y

    def temporal_discre_3rdorder_adamsbashforth(self, convective_prev, convective_prev2):
        convective_x = 23/12 * self.convective_x - 16/12 * convective_prev.convective_x + 5/12 * convective_prev2.convective_x
        convective_y = 23/12 * self.convective_y - 16/12 * convective_prev.convective_y + 5/12 * convective_prev2.convective_y
        return convective_x, convective_y


