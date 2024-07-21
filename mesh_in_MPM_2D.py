# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 16:20:03 2022

@author: user
"""

import numpy as np

class Material(object):
    def __init__(self, rho, K, G):
        self.rho = rho
        self.K = K
        self.G = G
        self.E = 9.0*K*G/(3.0*K+G)
        self.v = 0.5*(3.0*K-2.0*G)/(3.0*K+G) # poisson

        # plane stress condition
        self.E1 = self.E/(1-self.v*self.v)
        self.E2 = self.v*self.E/(1-self.v*self.v)
        self.D = np.array([
                            [self.E1, self.E2,    0.0],
                            [self.E2, self.E1,    0.0],
                            [    0.0,     0.0, self.G]
                          ])

    def verify_time_step(self, Unit_Grid, dt):
        # P wave velocity -> vel_p
        self.M = self.K+4.0*self.G/3.0
        self.vel_p = np.sqrt(self.M/self.rho)
        self.dt_critical = Unit_Grid/(10*self.vel_p)
        print(f"MPM critical time step = {self.dt_critical}")
        if (dt > self.dt_critical):
            from sys import exit
            print("dt =", dt)
            print("dt_critical =", self.dt_critical)
            print("need reset MPM time step")
            exit()


class Element(object):
    def __init__(self, Unit_Grid):
        self.eid = 0
        self.n1  = None
        self.n2  = None
        self.n3  = None
        self.n4  = None
        self.L   = np.array([Unit_Grid, Unit_Grid])
        self.particles = [0]  

class Particle(object):
    def __init__(self, material, position):
        self.rho      = material.rho
        self.E        = material.E
        self.position = position
        
        self.pid              = 0
        self.Volume           = 0.0
        self.mass             = 0.0
        self.velocity         = np.array([0.0, 0.0])
        self.momentum         = np.array([0.0, 0.0])
        self.body_force       = np.array([0.0, 0.0])
        self.strain           = np.array([0.0, 0.0, 0.0])
        self.stress           = np.array([0.0, 0.0, 0.0]) # sigmx, sigmy, tauxy
        self.specific_stress  = np.array([0.0, 0.0, 0.0])
        self.strain_increment = np.array([0.0, 0.0, 0.0])

        self.Ni  = np.array([0.0, 0.0, 0.0, 0.0]) # isoparametric 4 nodes element
        self.dNi = np.array([[0.0, 0.0], 
                             [0.0, 0.0],
                             [0.0, 0.0],
                             [0.0, 0.0]]) # np.zeros((4, 2))
        self.element = None
        

class Node(object):
    def __init__ (self):
        self.nid          = 0
        self.position     = np.array([0.0, 0.0])
        self.mass         = 0.0
        self.velocity     = np.array([0.0, 0.0])
        self.momentum     = np.array([0.0, 0.0])
        self.acceleration = np.array([0.0, 0.0])
        self.f_int        = np.array([0.0, 0.0])
        self.f_ext        = np.array([0.0, 0.0])
        self.f_bc         = np.array([0.0, 0.0])
        self.f_c          = np.array([0.0, 0.0])
        self.normal       = np.array([0.0, 0.0])
        self.backup_normal = np.array([0.0, 0.0])
        self.centroid_position = np.array([0.0, 0.0])
        self.IsDetectFEM  = False


class Mesh(object):
    def __init__(self, vol, xp, xn, nelem, Unit_Grid, tot_elem_num, tot_node_num, tot_ptc_num, elastic):
        self.elastic      = elastic
        self.tot_elem_num = tot_elem_num
        self.tot_node_num = tot_node_num
        self.tot_ptc_num  = tot_ptc_num
        self.elements     = [0]*tot_elem_num
        self.nodes        = [0]*tot_node_num
        self.particles    = [0]*tot_ptc_num
        
        for p in range(tot_ptc_num):
            self.particles[p] = Particle(elastic, xp[p])
            self.particles[p].pid = p
            self.particles[p].Volume = vol[p]

        for i in range(tot_node_num):
            self.nodes[i]          = Node()
            self.nodes[i].nid      = i
            self.nodes[i].position = xn[i]
      
        # put nodes in element
        particles_id = np.zeros(tot_ptc_num)
        i = 0
        for e in range(tot_elem_num):
            if ((e % nelem) == 0 and (e != 0)):
                i += 1
            self.elements[e] = Element(Unit_Grid)
            self.elements[e].eid = e
            self.elements[e].n1 = self.nodes[e+i]
            self.elements[e].n2 = self.nodes[e+1+i]
            self.elements[e].n3 = self.nodes[e+2+nelem+i] # 1
            self.elements[e].n4 = self.nodes[e+1+nelem+i] # 2
            
            # put particles in element
            self.elements[e].particles = []
            ptc_in_elem_num = 0
            j = 0
            for ip in self.particles:
                Xp = ip.position[0]
                Yp = ip.position[1] 
                if ((self.elements[e].n1.position[0] <= Xp) and (Xp < self.elements[e].n2.position[0])):
                    if ((self.elements[e].n2.position[1] <= Yp) and (Yp < self.elements[e].n3.position[1])):
                        if (particles_id[j] == 0):
                            self.elements[e].particles.append(ip)
                            self.elements[e].particles[ptc_in_elem_num].element = self.elements[e]
                            ptc_in_elem_num += 1
                            particles_id[j] += 1
                j += 1
                        
    
    def imposeMPM_IC(self, vpo, epo, bpo):
        i = 0
        for ip in self.particles:
            ip.velocity[0] = vpo[i, 0]
            ip.velocity[1] = vpo[i, 1]
            
            ip.strain[0] = epo[i, 0]
            ip.strain[1] = epo[i, 1]
            ip.strain[2] = epo[i, 2]
            
            ip.body_force[0] = bpo[i, 0]
            ip.body_force[1] = bpo[i, 1]
            i += 1
            
    # def impose_IC(self, bc, vpo, epo, bpo, xp):
    #     # initial function
    #     if (epo[0] == "function"):
    #         rho = self.elastic.rho
    #         E = self.elastic.E
    #         E1 = self.elastic.E1
    #         E2 = self.elastic.E2
    #         h = self.elements[0].L[0]
    #         gx, gy = bpo
    #         g = -pow(gx*gx+gy*gy, 0.5)

    #         n1_1 = np.array([xp[0, 0], -1.5])
    #         n1_2 = np.array([xp[4, 0], -1.5])
    #         n1_3 = np.array([xp[8, 0], -1.5])
    #         n1_4 = np.array([xp[12, 0], -1.5])
            
    #         exl = np.array([1.0, 0.0])
    #         eyl = np.array([0.0, 1.0])
    #         p1 = bc.p1
    #         p2 = bc.p2
            
    #         det = 1/np.sqrt(pow(p2[0]-p1[0], 2) + pow(p2[1]-p1[1], 2))
    #         ex = det*np.array([p2[0]-p1[0], p2[1]-p1[1]])
    #         ey = det*np.array([-p2[1]+p1[1], p2[0]-p1[0]])
                
    #         T = np.array([[ex.dot(exl), ey.dot(exl)],
    #                       [ex.dot(eyl), ey.dot(eyl)]])
            
    #         new_n1_1 = T.dot(n1_1-p1)+p1
    #         new_n1_2 = T.dot(n1_2-p1)+p1
    #         new_n1_3 = T.dot(n1_3-p1)+p1
    #         new_n1_4 = T.dot(n1_4-p1)+p1
            
    #         xbc = np.sort(np.array([new_n1_4[0], new_n1_3[0], new_n1_2[0], new_n1_1[0]]))
    #         ybc = np.sort(np.array([new_n1_4[1], new_n1_3[1], new_n1_2[1], new_n1_1[1]]))
            
    #         def ey(p_x, x_bc, g):
    #             px, py = p_x
    #             xbc, ybc = x_bc
    #             dy = pow(pow(px-xbc, 2)+pow(py-ybc, 2), 0.5)
    #             return rho*g*(h-dy)/E
        
    #     i, k = 0, 0
    #     for ip in self.particles:
    #         if (epo[0] == "function"):
    #             p_x = ip.position
    #             e_y = ey(p_x, [xbc[i], ybc[i]], g)
    #             e_x = -E2*e_y/E1

    #         for d in range(self.dim):
    #             ip.velocity[d] = vpo[d]
    #             ip.body_force[d] = bpo[d]
            
    #         for j in range(3):
    #             if (epo[0] != "function"):
    #                 ip.strain[j] = epo[j]
                        
    #             if (epo[0] == "function"):
    #                 if (j == 0):
    #                     ej = e_x
    #                 elif (j == 1):
    #                     ej = e_y
    #                 else:
    #                     ej = 0.0
    #                 ip.strain[j] = ej
        
    #         if ((k+1) % 4 == 0):
    #             i += 1
    #         k += 1
    
    # def imposeFEM_IC(self, bc, vpo, epo, bpo, xp, theta):
    #     # initial function
    #     if (epo[0] == "function"):
    #         rho = self.elastic.rho
    #         E = self.elastic.E
    #         E1 = self.elastic.E1
    #         E2 = self.elastic.E2
    #         h = self.elements[0].L[0]
    #         gx, gy = bpo
    #         g = -pow(gx*gx+gy*gy, 0.5)

    #         n1_1 = np.array([xp[0, 0], -1.5])
    #         n1_2 = np.array([xp[4, 0], -1.5])
    #         n1_3 = np.array([xp[8, 0], -1.5])
    #         n1_4 = np.array([xp[12, 0], -1.5])
            
    #         from creat_object_in_MPM_2D import Transform_cood_FEM
    #         p1, p2, T = Transform_cood_FEM(xp, bc, thickness=0.25, theta=theta)
            
    #         new_n1_1 = T.dot(n1_1-p1)+p1
    #         new_n1_2 = T.dot(n1_2-p1)+p1
    #         new_n1_3 = T.dot(n1_3-p1)+p1
    #         new_n1_4 = T.dot(n1_4-p1)+p1
            
    #         xbc = np.sort(np.array([new_n1_4[0], new_n1_3[0], new_n1_2[0], new_n1_1[0]]))
    #         ybc = np.sort(np.array([new_n1_4[1], new_n1_3[1], new_n1_2[1], new_n1_1[1]]))
            
    #         def ey(p_x, x_bc, g):
    #             px, py = p_x
    #             xbc, ybc = x_bc
    #             dy = pow(pow(px-xbc, 2)+pow(py-ybc, 2), 0.5)
    #             return rho*g*(h-dy)/E
        
    #     i, k = 0, 0
    #     for ip in self.particles:
    #         if (epo[0] == "function"):
    #             p_x = ip.position
    #             e_y = ey(p_x, [xbc[i], ybc[i]], g)
    #             e_x = -E2*e_y/E1

    #         for d in range(self.dim):
    #             ip.velocity[d] = vpo[d]
    #             ip.body_force[d] = bpo[d]
            
    #         for j in range(3):
    #             if (epo[0] != "function"):
    #                 ip.strain[j] = epo[j]
                        
    #             if (epo[0] == "function"):
    #                 if (j == 0):
    #                     ej = e_x
    #                 elif (j == 1):
    #                     ej = e_y
    #                 else:
    #                     ej = 0.0
    #                 ip.strain[j] = ej
        
    #         if ((k+1) % 4 == 0):
    #             i += 1
    #         k += 1
  
class BC(object):
    def __init__(self, form, p1, p2):
        self.form = form
        # self.Unit_Grid = Unit_Grid
        self.p1 = p1
        self.p2 = p2
        self.vbc = np.array([0.0, 0.0])
        self.mu = 0.0
        
        dy = self.p2[1]-self.p1[1]
        dx = self.p2[0]-self.p1[0]
        if self.p2[0]-self.p1[0] != 0:
            self.theta = np.arctan(dy/dx)
            if (dx < 0) and (dy > 0):
                self.theta = np.pi+self.theta
            if (dx < 0) and (dy < 0):
                self.theta = np.pi+self.theta
            if (dx > 0) and (dy < 0):
                self.theta = 2.0*np.pi+self.theta
            if (dx < 0) and (dy == 0):
                self.theta = np.pi
        if self.p2[0]-self.p1[0] == 0:
            if dy > 0:
                self.theta = np.pi/2.0
            else:
                self.theta = 3.0*np.pi/2.0
        
    def set_normal(self, nbc):
        self.nbc = nbc
        
    def set_friction_coefficient(self, mu):
        self.mu = mu
        
    def set_function(self, *arg):
        f = lambda x : arg[0]*x*x+arg[1]*x+arg[1]
        return f
        
    def f(self, x):
        # a = self.nbc[0]
        # b = self.nbc[1]
        # if (b == 0.0):
        #     return self.p1[0]
        # else:
        #     c = -b*self.p1[1]-a*self.p1[0]
        #     return (1.0/b)*(-a*x-c)
        return 0.0
        
    def df(self, x):
        # a = self.nbc[0]
        # b = self.nbc[1]
        # if (b == 0.0):
        #     return 0.0
        # else:
        #     return (1.0/b)*(-a)
        return 0.0
    
    def ddf(self, x):
        return 0.0
    
    def re_transform_matrix(self, xo):
        x = (xo[0]-self.p1[0])*np.cos(self.theta)+(xo[1]-self.p1[1])*np.sin(self.theta)
        y = -(xo[0]-self.p1[0])*np.sin(self.theta)+(xo[1]-self.p1[1])*np.cos(self.theta)
        return (x, y)
    
    def transform_matrix(self, xl):
        x = (xl[0])*np.cos(self.theta)-(xl[1])*np.sin(self.theta)
        y = (xl[0])*np.sin(self.theta)+(xl[1])*np.cos(self.theta)
        return (x+self.p1[0], y+self.p1[1])
  
    def newton_method(self, x_ik):
        tol = 1e-10
        xo = x_ik[0]
        yo = x_ik[1]
        x = xo
        y = round(self.f(x), 10)

        num = 0
        while (True):
            Tt = np.array([1, self.df(x)])

            if ((yo-y) == 0.0):
                n = np.array([self.nbc[0], self.nbc[1]])
                if (x < xo):
                    n *= (-1.0)
                A, dA, B, dB = 0.0, 0.0, 0.0, 0.0 
            else:
                n = np.array([xo-x, yo-y])/np.sqrt((xo-x)*(xo-x)+(yo-y)*(yo-y))
                A = 1/np.sqrt((xo-x)*(xo-x)+(yo-y)*(yo-y))
                dA = -0.5*pow((xo-x)*(xo-x)+(yo-y)*(yo-y), -1.5)*(-2.0*(xo-x)-2.0*(yo-y)*self.df((x)))
                B = Tt[0]*(xo-x)+Tt[1]*(yo-y)
                dB = -Tt[0]+self.ddf(x)*(yo-y)+Tt[1]*(-self.df(x))

            N = Tt.dot(n)
            if (abs(N) <= tol):
                return x, y
 
            if (num > 100):
                return x, y
            
            dN = dA*B+A*dB
            if (dN != 0.0):
                x -= N/dN
            y = round(self.f(x), 10)
            
            num += 1
            
            


if __name__ == "__main__" :
    
    rho = 960.0
    K = 0.003E+9
    G = 0.0006E+9
    
    Unit_Grid = 0.25
    dt = 1e-4
    
    elastic = Material(rho, K, G)
    elastic.verify_time_step(Unit_Grid, dt)
    print("wave speed =", elastic.vel_p)
    print("dt critical =", elastic.dt_critical)
    
    bc1 = BC("slip", p1=np.array([1.5,  -2.5]), p2=np.array([1.5,  2.5]))
    bc2 = BC("slip", p1=np.array([-1.5, -2.5]), p2=np.array([-1.5, 2.5]))
    bc3 = BC("slip", p1=np.array([-2.5,  1.5]), p2=np.array([2.5,  1.5]))
    bc4 = BC("slip", p1=np.array([-2.5, -1.5]), p2=np.array([2.5*pow(3, 0.5)-2.5, 1.0]))
    bc1.set_normal(nbc=np.array([-1.0, 0.0]))
    bc2.set_normal(nbc=np.array([1.0,  0.0]))
    bc3.set_normal(nbc=np.array([0.0, -1.0]))
    bc4.set_normal(nbc=np.array([-1/2, pow(3, 0.5)/2])) # 
    bc4.set_friction_coefficient(mu=0.70)
    bc_array=[bc4] # bc1, bc2, bc3, 
    
    x_ik = np.array([-2.5, -1.5])
    import matplotlib.pyplot as plt
    plt.figure(figsize=(16, 9))
    for bc in bc_array:
        x_ikl = bc.re_transform_matrix(x_ik)
        x_bc = bc.transform_matrix(bc.newton_method(x_ikl))
        D = (x_ik-x_bc).dot(bc.nbc)

        if (D <= 0.25):
            print(bc.nbc, D, x_ik, x_bc)
        plt.plot([bc.p1[0], bc.p2[0]], [bc.p1[1], bc.p2[1]], "-", lw=5, color="tab:blue")
        plt.plot(x_bc[0], x_bc[1], "o", markersize=15, color="tab:red")
        plt.plot(x_ik[0], x_ik[1], "o", markersize=15, color="black")
    plt.grid(True)
    plt.axos('equal')
    plt.show()
    
    
        
    