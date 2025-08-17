import os
import config
import my_io
import particle_tracking
import smac
import numpy as np

DT = config.DT
NX = config.NX
NY = config.NY

def main():
    output_folder_particle = "particle_tracking"
    os.makedirs(output_folder_particle, exist_ok=True)

    t = 0
    data_prev  = smac.StepData()
    data_prev  = smac.set_boundary_condition(data_prev)
    my_io.write_vti_flow(data_prev)
    data_prev2 = None
    data_prev3 = None
    particles_new = particle_tracking.scatter_particles() 
    
    while t < config.T_END:
        t += DT
        data_current, data_prev, data_prev2, data_prev3 = smac.calc_onestep_implicit(data_prev, data_prev2, data_prev3)
        
        my_io.write_vti_flow(data_current)

        # 残差計算, 収束判定
        numerator   = 0
        denominator = 0
        eps = 1
        for i in range (NX):
            for j in range (NY):
                numerator += (data_prev.u[i][j] - data_prev2.u[i][j])**2 + (data_prev.v[i][j] - data_prev2.v[i][j])**2
                denominator += (data_prev2.u[i][j])**2 + (data_prev2.v[i][j])**2
        if t > DT:
            eps = numerator/denominator
            if eps <= config.EPS_C:
                print("The flow field has reached steady state.")
                break

        U_average = 0
        countor   = 0
        for i in range(NX):
            for j in range(NY):
                U_average += np.sqrt(data_current.u[i][j]**2 + data_current.v[i][j]**2)
                countor   += 1
        U_average /= countor

        print(f"t = {t:.4f}", f"U_average = {U_average}", f"U_residual = {eps}")

        # 粒子trace
        particles = particles_new
        particles_new = []
        for particle in particles:
            particle_new = particle_tracking.rk4_method(particle, data_current)
            particles_new.append(particle_new)
        vti_filename_particle = os.path.join(output_folder_particle, f"t_{int(t/DT):04d}.vtk")
        my_io.write_vti_particle(vti_filename_particle, particles_new)

if __name__ == "__main__":
    main()
