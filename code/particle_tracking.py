#         o --- o --- o --- o   
#         |     |     |     |
#         |     |     |     |
#         o --- o --- o --- o
#         |     |     |     |
#         |  p  |     |     |
#         o --- o --- o --- o  
#         |     |     |     |
#         |     |     |     |
#         o --- o --- o --- o
#
# 
#        NX = NY = 4
#        DX = LX / (NX-1)
#        position of particle (p)
#        x_int = 0
#        y_int = 1

import random
import numpy as np
import config

DX = config.DX
DY = config.DY
LX = config.LX
LY = config.LY
DT = config.DT
NUM_OF_PARTICLE = config.NUM_OF_PARTICLE

class Particle:
    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y
        self.x_int  = int( x / DX )    # particleが所属するcellの, 左下の格子点のi番号, j番号
        self.y_int  = int( y / DY )    # 0 ~ NX-2 の NX-1 個の数になりうる
        # normalize at the cell
        while int (x / DX) > 0 :
            x -= DX
        while int (y / DY) > 0 :
            y -= DY
        self.x_norm = x / DX
        self.y_norm = y / DY

    def calc_velocity_at_the_cell(self, u, v):
        u_sw = u[self.x_int][self.y_int]
        u_se = u[self.x_int + 1][self.y_int]
        u_nw = u[self.x_int][self.y_int + 1]
        u_ne = u[self.x_int + 1][self.y_int + 1]
        v_sw = v[self.x_int][self.y_int]
        v_se = v[self.x_int + 1][self.y_int]
        v_nw = v[self.x_int][self.y_int + 1]
        v_ne = v[self.x_int + 1][self.y_int + 1]

        self.u = self.x_norm * self.y_norm * u_ne + ( 1 - self.x_norm) * self.y_norm * u_nw + self.x_norm * ( 1 - self.y_norm) * u_se + (1 - self.x_norm) * (1 - self.y_norm) * u_sw
        self.v = self.x_norm * self.y_norm * v_ne + ( 1 - self.x_norm) * self.y_norm * v_nw + self.x_norm * ( 1 - self.y_norm) * v_se + (1 - self.x_norm) * (1 - self.y_norm) * v_sw

def scatter_particles():
    particles_x = [random.uniform(0, LX) for _ in range(NUM_OF_PARTICLE)]
    particles_y = [random.uniform(0, LY) for _ in range(NUM_OF_PARTICLE)]
    particles  = []
    for i in range(len(particles_x)):
        particle = Particle(i, particles_x[i], particles_y[i])
        particles.append(particle)
    return particles

def boundary_condition(particle):
    if particle.x < 0:
        particle.x = 0       
    if particle.x >= LX:
        particle.x = LX - 0.00001*DX
    if particle.y < 0:
        particle.y = 0    
    if particle.y >= LY:
        particle.y = LY - 0.00001*DY
    particle_new = Particle(particle.id, particle.x, particle.y)
    return particle_new

def rk4_method(particle, data):
    particle.calc_velocity_at_the_cell(data.u, data.v)

    particle_tmp1 = Particle(particle.id, particle.x + DT*particle.u/2, particle.y + DT*particle.v/2)
    particle_tmp1 = boundary_condition(particle_tmp1)
    particle_tmp1.calc_velocity_at_the_cell(data.u, data.v)

    particle_tmp2 = Particle(particle.id, particle.x + DT*particle_tmp1.u/2, particle.y + DT*particle_tmp1.v/2)
    particle_tmp2 = boundary_condition(particle_tmp2)
    particle_tmp2.calc_velocity_at_the_cell(data.u, data.v)

    particle_tmp3 = Particle(particle.id, particle.x + DT*particle_tmp2.u, particle.y + DT*particle_tmp2.v)
    particle_tmp3 = boundary_condition(particle_tmp3)
    particle_tmp3.calc_velocity_at_the_cell(data.u, data.v)

    x_new = particle.x + DT*particle.u / 6 + DT * particle_tmp1.u/3 + DT*particle_tmp2.u / 3 + DT*particle_tmp3.u/6
    y_new = particle.y + DT*particle.v / 6 + DT * particle_tmp1.v/3 + DT*particle_tmp2.v / 3 + DT*particle_tmp3.v/6 

    particle_new = Particle(particle.id, x_new, y_new)
    particle_new = boundary_condition(particle_new)

    return particle_new