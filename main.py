import os
import flow
import config
import my_io
import particle_tracking
import numpy as np

DT = config.DT
NX = config.NX
NY = config.NY

def main():
    output_folder_flow = "flowfield"
    os.makedirs(output_folder_flow, exist_ok=True)
    output_folder_particle = "particle_tracking"
    os.makedirs(output_folder_particle, exist_ok=True)

    t    = 0
    u    = np.zeros((NX, NY))
    v    = np.zeros((NX, NY))
    p    = np.zeros((NX, NY)) 

    particles_new = particle_tracking.scatter_particles() 

    while t < config.T_END:
        u, v, p = flow.calc_onestep(u, v, p)
        t += DT
        print(f"t={t:.4f}")
        vti_filename_flow = os.path.join(output_folder_flow, f"t_{int(t/DT):04d}.vti")
        my_io.write_vti_flow(vti_filename_flow, u, v, p )

        #TODO: U magnitudeの可視化を含める
        #TODO: ここに残差計算コンソール表示を追加する
        
        particles = particles_new
        particles_new = []
        for particle in particles:
            particle_new = particle_tracking.rk4_method(particle, u, v)
            particles_new.append(particle_new)
        vti_filename_particle = os.path.join(output_folder_particle, f"t_{int(t/DT):04d}.vtk")
        my_io.write_vti_particle(vti_filename_particle, particles_new)

if __name__ == "__main__":
    main()
