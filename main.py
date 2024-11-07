import os
import time
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util import numpy_support

output_dir = "output_cavity"
os.makedirs(output_dir, exist_ok=True)

Lx = 1 
Ly = 1 
U = 1       # velocity of ceiling wall
Re = 1      #  Reynolds number 
Nx = 32     # grid counts in x axis
Ny = 32     #                y
t_end = 1.0 # total time
dt = 0.0001 #delta t
eps_c = 0.00000001  # Threshold

Ncells = Nx * Ny 
dx = Lx/Nx
dy = Ly/Ny
nu = U*Lx/Re 

def write_vti_file(filename, u_visu, v_visu, p_visu, KE_visu, Nx, Ny):

    imageData = vtk.vtkImageData()
    imageData.SetDimensions(Nx, Ny, 1)  

    u_vtk = numpy_support.numpy_to_vtk(u_visu.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    v_vtk = numpy_support.numpy_to_vtk(v_visu.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    p_vtk = numpy_support.numpy_to_vtk(p_visu.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)
    KE_vtk = numpy_support.numpy_to_vtk(KE_visu.ravel(), deep=True, array_type=vtk.VTK_DOUBLE)

    imageData.GetPointData().AddArray(u_vtk)
    imageData.GetPointData().AddArray(v_vtk)
    imageData.GetPointData().AddArray(p_vtk)
    imageData.GetPointData().AddArray(KE_vtk)


    u_vtk.SetName("u")
    v_vtk.SetName("v")
    p_vtk.SetName("p")
    KE_vtk.SetName("KE")

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(imageData)
    writer.Write()

u = np.zeros((Nx+1,Ny))
up = np.zeros((Nx+1,Ny))   # u prediction  (SMAC method)
uold =np.zeros((Nx+1,Ny))   # for convergence check  
v = np.zeros((Nx, Ny+1))
vp = np.zeros((Nx, Ny+1))     # v prediction (SMAC method)
vold =np.zeros((Nx,Ny+1))    # for convergence check
phi = np.zeros(Nx*Ny)          #pressure modification
M = np.zeros((Ncells, Ncells))    #Coefficient matrix of the pressure Poisson equation
b = np.zeros(Ncells) 
A1=np.zeros((Nx-1,Ny)) 
B1=np.zeros((Nx-1,Ny))
A2=np.zeros((Nx,Ny-1))
B2=np.zeros((Nx,Ny-1))
p=np.zeros((Nx,Ny))    
dpdx=np.zeros((Nx-1,Ny))
dpdy=np.zeros((Nx,Ny-1))
uxx=np.zeros((Nx-1,Ny))
uyy=np.zeros((Nx-1,Ny))
vxx=np.zeros((Nx,Ny-1))
vyy=np.zeros((Nx,Ny-1))
KE=np.zeros((Nx,Ny))

#prepare M
for i in range (Nx):
    for j in range (Ny):
        PP = i + j * Nx
        SS = i + (j - 1) * Nx
        WW = (i - 1) + j * Nx
        EE = (i + 1) + j * Nx
        NN = i + (j + 1) * Nx
        if j >=1:
            M[PP][SS] = 1 / (dy * dy)
            M[PP][PP] -= 1 / (dy * dy)
        if i >= 1:
            M[PP][WW] = 1 / (dx * dx)
            M[PP][PP] -= 1 / (dx * dx)
        if i <= Nx - 2:
            M[PP][EE] = 1 / (dx * dx)
            M[PP][PP] -= 1 / (dx * dx)
        if j <= Ny - 2:
            M[PP][NN] = 1 / (dy * dy)
            M[PP][PP] -= 1 / (dy * dy)

t = dt
# main loop
while t <= t_end:
    for i in range (Nx+1):
        for j in range (Ny):
            uold[i][j] = u[i][j]
    for i in range (Nx):
        for j in range (Ny+1):
            vold[i][j] = v[i][j]

############  calculate up #########################
    # calculate A1
    for j in range (1, Ny-1):
        A1[0][j] = -(1/dx)*(-(u[0][j]/2)**2+((u[0][j]+u[1][j])/2)**2) \
                    -(1/dy)*(-(v[0][j-1]+v[1][j-1])*(u[0][j-1]+u[0][j])/4+(v[0][j]+v[1][j])*(u[0][j+1]+u[0][j])/4)
        A1[Nx-2][j] = -(1/dx)*(-((u[Nx-2][j]+u[Nx-3][j])/2)**2+(u[Nx-2][j]/2)**2) \
                    -(1/dy)*(-(v[Nx-1][j-1]+v[Nx-2][j-1])*(u[Nx-2][j-1]+u[Nx-2][j])/4+(v[Nx-1][j]+v[Nx-2][j])*(u[Nx-2][j+1]+u[Nx-2][j])/4)

    for i in range (1,Nx-2):                                           
        A1[i][0] = -(1/dx)*(-((u[i-1][0]+u[i][0])/2)**2+((u[i+1][0]+u[i][0])/2)**2) \
                    -(1/dy)*((v[i][0]+v[i+1][0])*(u[i][1]+u[i][0])/4)
        A1[i][Ny-1] = -(1/dx)*(-((u[i][Ny-1]+u[i-1][Ny-1])/2)**2+((u[i+1][Ny-1]+u[i][Ny-1])/2)**2) \
                    -(1/dy)*(-(v[i][Ny-2]+v[i+1][Ny-2])*(u[i][Ny-2]+u[i][Ny-1])/4)             

    A1[0][0] = -(1/dx)*(((-u[0][0]/2)**2)+((u[1][0]+u[0][0])/2)**2) \
                - (1/dy)*(v[1][0]+v[0][0])*(u[0][1]+u[0][0])/4                        
    A1[0][Ny-1] = -(1/dx)*(-((u[0][Ny-1])/2)**2+((u[1][Ny-1]+u[0][Ny-1])/2)**2) \
                - (1/dy)*(-(v[0][Ny-2]+v[1][Ny-2])*(u[0][Ny-2]+u[0][Ny-1])/4)        
    A1[Nx-2][0] = -(1/dx)*(-((u[Nx-3][0]+u[Nx-2][0])/2)**2+((u[Nx-2][0])/2)**2) \
                -(1/dy)*((u[Nx-2][1]+u[Nx-2][0])*(v[Nx-2][0]+v[Nx-1][0])/4)         
    A1[Nx-2][Ny-1] = -(1/dx)*(-((u[Nx-3][Ny-1]+u[Nx-2][Ny-1])/2)**2+((u[Nx-2][Ny-1])/2)**2) \
                        -(1/dy)*(-(v[Nx-2][Ny-2]+v[Nx-1][Ny-2])*(u[Nx-2][Ny-2]+u[Nx-2][Ny-1])/4)

    for i in range (1,Nx-2):
        for j in range (1, Ny-1):
            A1[i][j] = -(1/dx)*(-((u[i-1][j]+u[i][j])/2)**2+((u[i][j]+u[i+1][j])/2)**2) \
                        - (1/dy)*(-(u[i][j-1]+u[i][j])*(v[i][j-1]+v[i+1][j-1])/4+(u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])/4)     
            
    # calculate dpdx
    for i in range (Nx-1):
        for j in range (Ny):
            dpdx[i][j] = (p[i+1][j]-p[i][j])/dx  

    ##### calculate B1
    # calculate u_xx
    for j in range (Ny):
        uxx[0][j]=(-2*u[0][j]+u[1][j])/((dx)**2)                  
        uxx[Nx-2][j]=(u[Nx-3][j]-2*u[Nx-2][j])/((dx)**2)      
    for i in range (1,Nx-2):
        for j in range (Ny):
            uxx[i][j] = (u[i-1][j]-2*u[i][j]+u[i+1][j])/(dx**2)         
    #calculate u_yy

    for i in range (Nx-1):
        uyy[i][0] = (-12*u[i][0]+4*u[i][1])/(3*((dy)**2))

        uyy[i][Ny-1] = (4*u[i][Ny-2]-12*u[i][Ny-1]+8*U)/(3*((dy)**2))    

    for i in range (Nx-1):
        for j in range (1,Ny-1):
            uyy[i][j] = (u[i][j-1]-2*u[i][j]+u[i][j+1])/((dy)**2) 

    B1=nu*(uxx+uyy)  
    ####

    
    # calculate up with 1st order Euler method
    for i in range (Nx-1):
        for j in range (Ny):
                up[i][j]=u[i][j]+dt*(A1[i][j]-dpdx[i][j]+B1[i][j])  
########################################

########## calculate vp  ##############################
    # calculate A2
    for j in range (1,Ny-2):
        A2[0][j] = -(1/dx)*((u[0][j]+u[0][j+1])*(v[1][j]+v[0][j])/4)   -   (1/dy)*(-((v[0][j-1]+v[0][j])/2)**2+((v[0][j+1]+v[0][j])/2)**2)  

    for j in range (1,Ny-2):
        A2[Nx-1][j] = -(1/dx)*(-(u[Nx-2][j]+u[Nx-2][j+1])*(v[Nx-1][j]+v[Nx-2][j])/4)  -  (1/dy)*(-((v[Nx-1][j-1]+v[Nx-1][j])/2)**2+((v[Nx-1][j+1]+v[Nx-1][j])/2)**2)  

    for i in range (1,Nx-1):
        A2[i][0] = -(1/dx)*(-(u[i-1][1]+u[i-1][0])*(v[i-1][0]+v[i][0])/4+(u[i][1]+u[i][0])*(v[i+1][0]+v[i][0])/4) \
                    -(1/dy)*(-((v[i][0])/2)**2+((v[i][1]+v[i][0])/2)**2)                               

    for i in range (1,Nx-1):
        A2[i][Ny-2] = -(1/dx)*(-(u[i-1][Ny-1]+u[i-1][Ny-2])*(v[i-1][Ny-2]+v[i][Ny-2])/4+(u[i][Ny-1]+u[i][Ny-2])*(v[i+1][Ny-2]+v[i][Ny-2])/4) \
                        -(1/dy)*(-((v[i][Ny-3]+v[i][Ny-2])/2)**2+((v[i][Ny-2])/2)**2)

    A2[0][0] = -(1/dx)*((u[0][0]+u[0][1])*(v[1][0]+v[0][0])/4)  -  (1/dy)*(-((v[0][0])/2)**2+((v[0][0]+v[0][1])/2)**2)          
    A2[0][Ny-2] = -(1/dx)*((u[0][Ny-1]+u[0][Ny-2])*(v[1][Ny-2]+v[0][Ny-2])/4)  -  (1/dy)*(-((v[0][Ny-2]+v[0][Ny-3])/2)**2+((v[0][Ny-2])/2)**2)      
    A2[Nx-1][0] = -(1/dx)*(-(u[Nx-2][0]+u[Nx-2][1])*(v[Nx-2][0]+v[Nx-1][0])/4) \
                -(1/dy)*(-((v[Nx-1][0])/2)**2+((v[Nx-1][1]+v[Nx-1][0])/2)**2)                       
    A2[Nx-1][Ny-2] = -(1/dx)*(-(u[Nx-2][Ny-1]+u[Nx-2][Ny-2])*(v[Nx-2][Ny-2]+v[Nx-1][Ny-2])/4) - (1/dy)*(-((v[Nx-1][Ny-3]+v[Nx-1][Ny-2])/2)**2+((v[Nx-1][Ny-2])/2)**2)
                                                                                                    
    for i in range (1,Nx-1):
        for j in range (1,Ny-2):
            A2[i][j] = -(1/dx)*(-(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[j][j])/4+(u[i][j+1]+u[i][j])*(v[i+1][j]+v[i][j])/4) \
                        - (1/dy)*(-((v[i][j-1]+v[i][j])/2)**2+((v[i][j]+v[i][j+1])/2)**2)                   
    #

    #calculate dpdy
    for i in range (Nx):
        for j in range (Ny-1):
            dpdy[i][j] = (p[i][j+1]-p[i][j])/dy    
                
    #### calculate B2
    # calcualte vxx
    for j in range (Ny-1):
        vxx[0][j] = (-12*v[0][j]+4*v[1][j])/(3*((dx)**2))
        vxx[Nx-1][j] = (4*v[Nx-2][j]-12*v[Nx-1][j])/(3*((dx)**2))
    for i in range (1,Nx-1):
        for j in range (Ny-1):
            vxx[i][j] = (v[i-1][j]-2*v[i][j]+v[i+1][j])/(dx**2)
    # calculate vyy
    for i in range (Nx):
        vyy[i][0] = (-2*v[i][0]+v[i][1])/(dy**2)
        vyy[i][Ny-2] = (v[i][Ny-3]-2*v[i][Ny-2])/(dy**2) 
    for i in range (Nx):
        for j in range (1,Ny-2):
            vyy[i][j] = (v[i][j-1]-2*v[i][j]+v[i][j+1])/(dy**2)

    for i in range (Nx):
        for j in range (Ny-1):
            B2[i][j] = nu*(vxx[i][j]+vyy[i][j])
    ####

    # 1st order Euler method
    for i in range (Nx):
        for j in range (Ny-1):
            vp[i][j] = v[i][j]+dt*(A2[i][j]-dpdy[i][j]+B2[i][j])
    ################


    # Right side of the pressure equation : b
    for i in range (1,Nx-1):
        for j in range (1,Ny-1):
            b[i+j*Nx]=((-up[i-1][j]+up[i][j])/dx+(-vp[i][j-1]+vp[i][j])/dy)/dt 

    for i in range (1,Nx-1):
        b[i] = ((-up[i-1][0]+up[i][0])/dx+vp[i][0]/dy)/dt
        b[i+(Ny-1)*Nx] = ((-up[i-1][Ny-1]+up[i][Ny-1])/dx-vp[i][Ny-2]/dy)/dt 
                                                                  
    for j in range (1,Ny-1):
        b[j*Nx] = (up[0][j]/dx+(-vp[0][j-1]+vp[0][j])/dy)/dt
        b[Nx-1+j*Nx] = (-up[Nx-2][j]/dx+(-vp[Nx-1][j-1]+vp[Nx-1][j])/dy)/dt 
                                                                  
    b[0] = (up[0][0]/dx+vp[0][0]/dy)/dt   

    b[(Ny-1)*Nx] = (up[0][Ny-1]/dx-vp[0][Ny-2]/dy)/dt       

    b[Nx-1] = (-up[Nx-2][0]/dx+vp[Nx-1][0]/dy)/dt          

    b[Nx-1+(Ny-1)*Nx] = (-up[Nx-2][Ny-1]/dx-vp[Nx-1][Ny-2]/dy)/dt  

    ################    matrix calculation      ##############################
    phi=np.linalg.solve(M,b)
    ############

    # modify velocity and pressure 
    for i in range (Nx):
        for j in range (Ny):
            p[i][j] += phi[i+j*Nx]                 

    for i in range (Nx-1):
        for j in range (Ny):
            u[i][j]= up[i][j]-dt*(phi[i+1+j*Nx]-phi[i+j*Nx])/dx    
                                                                    

    for i in range (Nx):
        for j in range(Ny-1):
            v[i][j] = vp[i][j] -dt*(phi[i+(j+1)*Nx]-phi[i+j*Nx])/dy 

    for i in range (Nx):
        for j in range (Ny):
            KE[i][j] = 0.5*((u[i][j])**2 + (v[i][j])**2)

    # output
    u_visu=np.zeros((Nx,Ny))
    v_visu=np.zeros((Nx,Ny))
    p_visu=np.zeros((Nx,Ny))
    KE_visu=np.zeros((Nx,Ny))

    for i in range (Nx):
        for j in range (Ny):
            u_visu[i][j] = u[j][i]              
            v_visu[i][j] = v[j][i]              
            p_visu[i][j] = p[j][i]
            KE_visu[i][j] = KE[j][i]

    vti_filename = os.path.join(output_dir, f"t_{int(t*10000):04d}.vti")
    write_vti_file(vti_filename, u_visu, v_visu, p_visu, KE_visu, Nx, Ny)

    print("t=", f"{t:.4f}", end='  ')
    print("B1[Nx-2][Ny-1]=",f"{B1[Nx-2][Ny-1]:.2f}")

    ######### for convergence check ################
    numerator =0
    denominator =0
    eps = eps_c+1

    for i in range (Nx):
        for j in range (Ny):
            numerator += (u[i][j]-uold[i][j])**2 + (v[i][j]-vold[i][j])**2
            denominator += (uold[i][j])**2 + (vold[i][j])**2

    if t > dt:
        eps = numerator/denominator
        if eps <= eps_c:
            print("The flow field has reached steady state.")
            break

    t += dt
