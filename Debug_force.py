# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 16:40:31 2022

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt

NUM = 75

def input_ptc_data(FILE_PATH, N, tot_ptc_num):
    ptc_data = [0]*N
    with open(FILE_PATH, 'r') as file:
        i, r = 0, 0
        j = 0
        for text in file:
            if (r != 0):
                if (i % tot_ptc_num == 0):
                    ptc_data[int(i/tot_ptc_num)] = np.zeros((tot_ptc_num, 14))
                    j = 0
                else:
                    j += 1
                    
                text = text.split()
                for k in range(13):
                    if (k <= 1):
                        f = int(text[k])
                    else:
                        f = float(text[k])
                    ptc_data[int(i/tot_ptc_num)][j, k] = f
                i += 1
            r += 1
    return ptc_data

def input_node_data(FILE_PATH):
    node_data = []
    with open(FILE_PATH, 'r') as ifile:
        i = 0
        for text in ifile:
            text = text.split()
            if (i != 0):
                text[0] = int(text[12])
                text[1] = int(text[13])
                for j in range(2, 17):
                    text[j] = float(text[j])
                node_data.append(text)       
            i += 1
    return node_data

def deinput_ptc_data(func):
    def factory(FILE_PATH, N, tot_ptc_num):
        ptc_data = [0]*N
        with open(FILE_PATH, 'r') as file:
            i, r = 0, 0
            j = 0
            for text in file:
                if (r != 0):
                    if (i % tot_ptc_num == 0):
                        ptc_data[int(i/tot_ptc_num)] = np.zeros((tot_ptc_num, 14))
                        j = 0
                    else:
                        j += 1
                        
                    text = text.split()
                    for k in range(13):
                        if (k <= 1):
                            f = int(text[k])
                        else:
                            f = float(text[k])
                        ptc_data[int(i/tot_ptc_num)][j, k] = f
                    i += 1
                r += 1
        return func(ptc_data, N, tot_ptc_num)
    return factory

def deinput_node_data(func):
    def factory(FILE_PATH, N):
        node_data = []
        with open(FILE_PATH, 'r') as ifile:
            i = 0
            for text in ifile:
                text = text.split()
                if (True):
                    pass
                    # text[0] = int(text[2])
                    # text[1] = int(text[3])
                    for j in range(14):
                        text[j] = float(text[j])
                    node_data.append(text)       
                i += 1
        return func(node_data, N)
    return factory

def plot_node_boundary_MPM(ax, fig_id, A_line, B_line, TOT_line):
    ax.plot(tns[::1], A_line[::4], "-o", lw=5, color="tab:blue", markersize=13)  
    # ax.plot(tns[::NUM], B_line, "-s", lw=5, color="tab:green", markersize=9) 
    # ax.plot(tns, TOT_line, "-s", lw=5, color="tab:orange", markersize=13)
    ax.plot(-10, A_line[0], "-o", lw=5, color="tab:blue", label="MPM obj A", markersize=20)
    # ax.plot(-10, B_line[0], "-s", lw=5, color="tab:green", label="MPM obj B", markersize=20)
    # ax.plot(-10.0, 0.0, "-s", lw=5, color="tab:orange", label="MPM total", markersize=13)
            
    if (fig_id == 0): 
        ax.set_ylabel("$f_{bct}$ (N)", fontsize=25)
        ax.set_title("boundary force", fontsize=25)
    if (fig_id == 1):
        ax.legend(fontsize=25, loc=1, bbox_to_anchor=(0.9, 1.5), ncol=3)
        ax.set_ylabel("$f_{bcn}$ (N)", fontsize=25)
        ax.set_xlabel("time (s)", fontsize=25)
    
    ax.set_xlim(tns[0], tns[N-1])
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    plt.rc('font', size=20)
    
    ax.grid(True)
    # plt.show()

    
def plot_particle_vel_MPM(ax, fig_id, A_line, B_line, TOT_line):
    ax.plot(tns[::NUM], np.round(A_line, 10), "-o", lw=5, color="tab:blue", markersize=13)  
    ax.plot(tns[::NUM], np.round(B_line, 10), "-s", lw=5, color="tab:green", markersize=13) 
    ax.plot(-10, A_line[0], "-o", lw=5, color="tab:blue", label="MPM obj A", markersize=20)
    ax.plot(-10, B_line[0], "-s", lw=5, color="tab:green", label="MPM obj B", markersize=20)
    
    ax.set_xlim(tns[0], tns[N-1])
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    plt.rc('font', size=20)
    
    ax.grid(True)
    
    if (fig_id == 2): 
        ax.set_ylabel("$v_{pt}$ (m/s)", fontsize=25)
        ax.set_title("particle velocity", fontsize=25)
    if (fig_id == 3):
        ax.legend(fontsize=25, loc=1, bbox_to_anchor=(0.8, 1.5))
        ax.set_ylabel("$v_{pn}$ (m/s)", fontsize=25)   
        ax.set_xlabel("time (s)", fontsize=25)

@deinput_ptc_data
def PerParticleVelovity(ptc_data, N, tot_ptc_num):
    vpx = np.zeros((N, tot_ptc_num))
    vpy = np.zeros((N, tot_ptc_num))
    for i, ndata in enumerate(ptc_data):
        for j, data in enumerate(ndata):
            vpx[i, j] = data[5]
            vpy[i, j] = data[6]
    return vpx, vpy

@deinput_ptc_data
def PerParticlePosition(ptc_data, N, tot_ptc_num):
    xp = np.zeros((N, tot_ptc_num))
    yp = np.zeros((N, tot_ptc_num))
    for i, ndata in enumerate(ptc_data):
        for j, data in enumerate(ndata):
            xp[i, j] = data[3]
            yp[i, j] = data[4]
    return xp, yp

@deinput_ptc_data
def center_of_mass_velovity(ptc_data, N, tot_ptc_num):
    v_cm = np.zeros((N, 2))
    n = 0
    for ndata in ptc_data:
        M = 0
        for data in ndata:
            mp = data[2]
            vpx = data[5]
            # if (n == 0):
            #     print(data, "\n")
            vpy = data[6]
            M += mp
            v_cm[n, 0] += mp*vpx
            v_cm[n, 1] += mp*vpy
            
        if (M != 0):
            v_cm[n, 0] /= M
            v_cm[n, 1] /= M
        n += 1
    return v_cm[:, 0], v_cm[:, 1]

@deinput_ptc_data
def vel_avg(ptc_data, N, tot_ptc_num):
    avg = np.zeros((N, 2))
    for i in range(N):
        vpx = 0
        vpy = 0
        for j in range(tot_ptc_num):
            vpx += ptc_data[i][j, 5]
            vpy += ptc_data[i][j, 6]
        avg[i, 0] = vpx/tot_ptc_num
        avg[i, 1] = vpy/tot_ptc_num
    return avg

def stress_avg(ptc_data):
    avg = np.zeros((N, 3))
    for i in range(N):
        sx = 0
        sy = 0
        sxy = 0
        for j in range(tot_ptc_num):
            sx += ptc_data[i][j, 7]
            sy += ptc_data[i][j, 8]
            sxy += ptc_data[i][j, 9]
        avg[i, 0] = sx/tot_ptc_num
        avg[i, 1] = sy/tot_ptc_num
        avg[i, 2] = sxy/tot_ptc_num
    return avg

def sum_mass(node_data, N):
    s_mass = np.zeros(N)
    i, j = 0, 0
    mi = 0
    for data in node_data:
        tns = data[0]
        if (tns == j):
            mi += data[2]
        if (tns != j):
            s_mass[j] = mi
                
            mi = 0
            mi += data[2]
            j += 1
        if (i == (len(node_data)-1)):
            s_mass[j] = mi
        i += 1
    return s_mass

@deinput_node_data
def f_avg(node_data, N):
    avg = np.zeros((N, 8))
    i = 0
    j = 0
    num = 0
    fintx = 0.0
    finty = 0.0
    fextx = 0.0
    fexty = 0.0
    fbcx = 0.0
    fbcy = 0.0
    fcx = 0.0
    fcy = 0.0
    for data in node_data:
        if (i != data[0]):
            avg[i, 0] += fintx #/num
            avg[i, 1] += finty #/num
            avg[i, 2] += fextx #/num
            avg[i, 3] += fexty #/num
            avg[i, 4] += fbcx #/num
            avg[i, 5] += fbcy #/num
            avg[i, 6] += fcx
            avg[i, 7] += fcy
            fintx = 0.0
            finty = 0.0
            fextx = 0.0
            fexty = 0.0
            fbcx = 0.0
            fbcy = 0.0
            fcx = 0.0
            fcy = 0.0
            i += 1 
            num = 0
        # fintx += data[9]
        # finty += data[10]
        # fextx += data[11]
        # fexty += data[12]
        fbcx += data[12]
        fbcy += data[13]
        # fcx += data[15]
        # fcy += data[16]
        if (j == (node_data.__len__()-1)):
            avg[i, 0] += fintx #/(num+1)
            avg[i, 1] += finty #/(num+1)
            avg[i, 2] += fextx #/(num+1)
            avg[i, 3] += fexty #/(num+1)
            avg[i, 4] += fbcx #/(num+1)
            avg[i, 5] += fbcy #/(num+1)
            avg[i, 6] += fcx
            avg[i, 7] += fcy
        j += 1
        num += 1   
    return avg

def plot_node_boundary(ax, fig_id, A_line, B_line, TOT_line):
    ax.plot(tns, A_line, "-", lw=5, color="tab:blue")  
    ax.plot(tns, B_line, "--", lw=5, color="tab:green") #
    ax.plot(tns, TOT_line, ":", lw=5, color="tab:orange")
    ax.plot(0.0, 0.0, "-", lw=5, color="tab:blue", label="obj A")
    ax.plot(0.0, 0.0, "--", lw=5, color="tab:green", label="obj B")
    ax.plot(0.0, 0.0, ":", lw=5, color="tab:orange", label="total")

    ax.set_xlim(tns[0], tns[N-1])
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    ax.yaxis.get_offset_text().set_fontsize(20)
    ax.legend(fontsize=25, loc=1, ncol=3)
    ax.grid(True)
    
    if (fig_id == 0): 
        ax.set_ylabel("$f_{bct}$ (N)", fontsize=25)
        ax.set_title("boundary force", fontsize=25)

    if (fig_id == 1):
        ax.set_ylabel("$f_{bcn}$ (N)", fontsize=25)
        ax.set_xlabel("time (s)", fontsize=25)
        plt.show()
            
def plot_particle_vel(ax, fig_id, A_line, B_line, TOT_line):
    ax.plot(tns, np.round(A_line, 8), "-", lw=5, color="tab:blue")  
    ax.plot(tns, np.round(B_line, 8), "--", lw=5, color="tab:green") 
    ax.plot(tns, np.round(TOT_line, 8), ":", lw=5, color="tab:orange")
    ax.plot(0.0, 0.0, "-", lw=5, color="tab:blue", label="obj A")
    ax.plot(0.0, 0.0, "--", lw=5, color="tab:green", label="obj B")
    ax.plot(0.0, 0.0, ":", lw=5, color="tab:orange", label="total")
    
    ax.legend(fontsize=25, loc=1, ncol=3)
    ax.set_xlim(tns[0], tns[N-1])       
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    ax.yaxis.get_offset_text().set_fontsize(20)
    
    ax.grid(True)
    
    if (fig_id == 2): 
        
        ax.set_ylabel("$v_{pt}$ (m/s)", fontsize=25)
        ax.set_title("particle velocity", fontsize=25)
        
    if (fig_id == 3):
        ax.set_ylabel("$v_{pn}$ (m/s)", fontsize=25)
        ax.set_xlabel("time (s)", fontsize=25)
        plt.show()
        
    if (fig_id == 6):
        ax.set_ylabel("$\sigma_x$ (Pa)", fontsize=25)
        ax.set_title("particle stress", fontsize=25)
        
    if (fig_id == 7):
        ax.set_ylabel("$\sigma_y$ (Pa)", fontsize=25)
    if (fig_id == 8): 
        ax.set_ylabel("$\u03C4_{xy}$ (Pa)", fontsize=25)
        ax.set_xlabel("time (s)", fontsize=25)
        plt.show()

    
def plot_collision_force(ax, fig_id, A_line, B_line, TOT_line):
    ax.plot(tns, A_line, "-", lw=5, color="tab:blue")  
    ax.plot(tns, B_line, "--", lw=5, color="tab:green")
    ax.plot(tns, TOT_line, ":", lw=5, color="tab:orange")
    ax.plot(0.0, 0.0, "-", lw=5, color="tab:blue", label="obj A")
    ax.plot(0.0, 0.0, "--", lw=5, color="tab:green", label="obj B")
    ax.plot(0.0, 0.0, ":", lw=5, color="tab:orange", label="total")
    
    ax.set_xlim(tns[0], tns[N-1])
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    ax.yaxis.get_offset_text().set_fontsize(20)
    ax.legend(fontsize=25, loc=1, ncol=3)
    ax.grid(True)
            
    if (fig_id == 4): 
        ax.set_ylabel("$f_{cx}$ (N)", fontsize=25)
        ax.set_title("collision force", fontsize=25)
        
    if (fig_id == 5):
        ax.set_ylabel("$f_{cy}$ (N)", fontsize=25)
        ax.set_xlabel("time (s)", fontsize=25)
        plt.show()

def plot_total_mass_MPM(ax, A_line, B_line):
    ax.plot(tns[::NUM], np.round(A_line, 10), "-o", lw=5, color="tab:blue", markersize=13)  
    ax.plot(tns[::NUM], np.round(B_line, 10), "-s", lw=5, color="tab:green", markersize=13) 
    ax.plot(-10, A_line[0], "-o", lw=5, color="tab:blue", label="MPM obj A", markersize=20)
    ax.plot(-10, B_line[0], "-s", lw=5, color="tab:green", label="MPM obj B", markersize=20)
    
    ax.set_xlim(tns[0], tns[N-1])
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    plt.rc('font', size=20)
    ax.grid(True)
    ax.set_ylabel("total mass (kg)", fontsize=25)
    ax.legend(fontsize=25, loc=1, ncol=2)
    ax.set_xlabel("time (s)", fontsize=25)
    
def plot_force(ax, fid, tns, f_int, f_ext, title):
    ax.plot(tns, f_int, "-", lw=5, color="tab:blue")  
    ax.plot(tns, f_ext, "--", lw=5, color="tab:green")
    # ax.plot(tns, f_int+f_ext, ":", lw=5, color="tab:orange")

    ax.plot(0.0, 0.0, "-", lw=5, color="tab:blue", label="internal")
    ax.plot(0.0, 0.0, "--", lw=5, color="tab:green", label="external")
    # ax.plot(0.0, 0.0, ":", lw=5, color="tab:orange", label="total")

    ax.set_xlim(tns[0], tns[N-1])
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    ax.yaxis.get_offset_text().set_fontsize(20)
    ax.legend(fontsize=25, loc=1, ncol=3)
    ax.grid(True)
    
    if (fid == 10):
        ax.set_title(title, fontsize=25)
        ax.set_ylabel("$f_x$ (N)", fontsize=25)
    
    if (fid == 11):
        ax.set_ylabel("$f_y$ (N)", fontsize=25)
        ax.set_xlabel("time (s)", fontsize=25)
        plt.show()


if __name__ == "__main__":

    from sys import exit
    
    func_case = 1
    alpha = 0.0
    theta = 0
    mu = 0.0

    FILE_PATH_A_N = r"C:\Users\lintim0622\source\repos\MPM\node_output.txt"
    # FILE_PATH_B_N = r"D:\MPM_Code\MPM_2D\f%dalpha%4.2ftheta%dmu%3.1fnode_Info_B.txt"%(func_case, alpha, theta, mu)

    FILE_PATH_A_P = r"C:\Users\lintim0622\source\repos\MPM\particle_output.txt"
    # FILE_PATH_B_P = r"D:\MPM_Code\MPM_2D\f%dalpha%4.2ftheta%dmu%3.1fparticle_Info_B.txt"%(func_case, alpha, theta, mu)

    EndTime = 0.5
    dt = 1e-4
    tot_node_num = 169 # , 1681
    tot_ptc_num = 16
    MPMRecStep = 1
    plot_form = "collision force" # nodal boundary force, particle velocity,  , particle stress, total force, sum mass
    
    if (theta == 0):
        theta = 0
    elif (theta == 30):
        theta = np.pi/6
    elif (theta == 45):
        theta = np.pi/4
    elif (theta == 90):
        theta = np.pi/2
    else:
        print("NO theta")

    if (plot_form == "particle velocity"):
        tns = np.arange(0.0, EndTime+dt*MPMRecStep/10.0, dt*MPMRecStep)
        N = int(round((EndTime/dt)/MPMRecStep+1, 8))
    else:
        tns = np.arange(0, EndTime+dt/10.0, dt)
        N = int(round(EndTime/dt, 7))+1
        
    plt.rcParams["font.family"] = "Times New Roman"
    if (plot_form == "nodal boundary force"):
        
        avg_fbcA = f_avg(FILE_PATH_A_N, 20004)[:, 4:6] #b f_avg(input_node_data(FILE_PATH_A_N), N)[:, 4:6]
        # avg_fbcB = f_avg(FILE_PATH_B_N, N)[:, 4:6] # f_avg(input_node_data(FILE_PATH_B_N), N)[:, 4:6]
        
        fbc_Ax, fbc_Ay = avg_fbcA[:, 0][::1], avg_fbcA[:, 1][::1]
        fbc_At = fbc_Ax*np.cos(theta) + fbc_Ay*np.sin(theta)
        fbc_An = -fbc_Ax*np.sin(theta) + fbc_Ay*np.cos(theta)
        
        # fbc_Bx, fbc_By = avg_fbcB[:, 0][::NUM], avg_fbcB[:, 1][::NUM]
        # fbc_Bt = fbc_Bx*np.cos(theta) + fbc_By*np.sin(theta)
        # fbc_Bn = -fbc_Bx*np.sin(theta) + fbc_By*np.cos(theta)
        
        fig, (ax1, ax2) = plt.subplots(2, 1)
        fig.subplots_adjust(hspace=0.7)
        # ax1.set_ylim(-1200.0, 1200.0)
        # ax2.set_ylim(300.0, 800.0)
        plot_node_boundary_MPM(ax1, 0, fbc_At, 0, fbc_At+0)
        plot_node_boundary_MPM(ax2, 1, fbc_An, 0, fbc_An+0)

    
    elif (plot_form == "particle velocity"):
        
        vpA_cmx, vpA_cmy = center_of_mass_velovity(FILE_PATH_A_P, N, tot_ptc_num) # center_of_mass_velovity(input_ptc_data(FILE_PATH_A_P, N, tot_ptc_num), N, tot_ptc_num)
        # vpB_cmx, vpB_cmy = center_of_mass_velovity(FILE_PATH_B_P, N, tot_ptc_num) # center_of_mass_velovity(input_ptc_data(FILE_PATH_B_P, N, tot_ptc_num), N, tot_ptc_num)

        fig, (ax1, ax2) = plt.subplots(2, 1)
        fig.subplots_adjust(hspace=0.5)
        
        vpA_t = vpA_cmy*np.sin(theta)+vpA_cmx*np.cos(theta)
        vpA_n = -vpA_cmy*np.cos(theta)+vpA_cmx*np.sin(theta)
        
        # vpB_t = vpB_cmy*np.sin(theta)+vpB_cmx*np.cos(theta)
        # vpB_n = vpB_cmy*np.cos(theta)-vpB_cmx*np.sin(theta)

        # plot_particle_vel(ax1, 2, vpA_t, vpB_t, vpA_t+vpB_t)
        # plot_particle_vel(ax2, 3, vpA_n, vpB_n, vpA_n+vpB_n)
        ax1.set_ylim(-6.0, 6.0)
        ax2.set_ylim(-6.0, 6.0)

    elif (plot_form == "particle stress"):
        
        avg_spA = stress_avg(input_ptc_data(FILE_PATH_A_P, N, tot_ptc_num))
        # avg_spB = stress_avg(input_ptc_data(FILE_PATH_B_P, N, tot_ptc_num))

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
        fig.subplots_adjust(hspace=0.5)
        # plot_particle_vel(ax1, 6, avg_spA[:, 0], avg_spB[:, 0], avg_spA[:, 0]+avg_spB[:, 0])
        # plot_particle_vel(ax2, 7, avg_spA[:, 1], avg_spB[:, 1], avg_spA[:, 1]+avg_spB[:, 1])
        # plot_particle_vel(ax3, 8, avg_spA[:, 2], avg_spB[:, 2], avg_spA[:, 2]+avg_spB[:, 2])

    elif (plot_form == "collision force"):
        
        avg_fcA = f_avg(input_node_data(FILE_PATH_A_N), N)[:, 6:]
        # avg_fcB = f_avg(input_node_data(FILE_PATH_B_N), N)[:, 6:]
        
        fig, (ax1, ax2) = plt.subplots(2, 1)
        fig.subplots_adjust(hspace=0.5)
        # plot_collision_force(ax1, 4, avg_fcA[:, 0], avg_fcB[:, 0], avg_fcA[:, 0]+avg_fcB[:, 0])
        # plot_collision_force(ax2, 5, avg_fcA[:, 1], avg_fcB[:, 1], avg_fcA[:, 1]+avg_fcB[:, 1])

    elif (plot_form == "total force"):
        
        avg_f_intA = f_avg(input_node_data(FILE_PATH_A_N), N)[:, :3]
        # avg_f_intB = f_avg(input_node_data(FILE_PATH_B_N), N)[:, :3]
        f_intAx = avg_f_intA[:, 0]
        f_intAy = avg_f_intA[:, 1]
        # f_intBx = avg_f_intB[:, 0]
        # f_intBy = avg_f_intB[:, 1]

        avg_f_extA = f_avg(input_node_data(FILE_PATH_A_N), N)[:, 2:4]
        # avg_f_extB = f_avg(input_node_data(FILE_PATH_B_N), N)[:, 2:4]
        f_extAx = avg_f_extA[:, 0]
        f_extAy = avg_f_extA[:, 1]
        # f_extBx = avg_f_extB[:, 0]
        # f_extBy = avg_f_extB[:, 1]
        
        fig, (ax1, ax2) = plt.subplots(2, 1)
        fig.subplots_adjust(hspace=0.5)
        plot_force(ax1, 10, tns, f_intAx, f_extAx, "obj A")
        plot_force(ax2, 11, tns, f_intAy, f_extAy, "obj A")
        
        fig, (ax1, ax2) = plt.subplots(2, 1)
        fig.subplots_adjust(hspace=0.5)
        # plot_force(ax1, 10, tns, f_intBx, f_extBx, "obj B")
        # plot_force(ax2, 11, tns, f_intBy, f_extBy, "obj B")
        
    elif (plot_form == "sum mass"):
        sum_MA = sum_mass(input_node_data(FILE_PATH_A_N), N)[::NUM]
        # sum_MB = sum_mass(input_node_data(FILE_PATH_B_N), N)[::NUM]
        
        fig, ax = plt.subplots()
        # plot_total_mass_MPM(ax, sum_MA, sum_MB)

        ax.set_ylim(10.0, 65.0)
        
    else:
        # A_data = input_ptc_data(FILE_PATH_A_P, N, tot_ptc_num)
        # v_pA = center_of_mass_velovity(A_data, N, tot_ptc_num)
        
        A_data = input_node_data(FILE_PATH_A_N)
        print("No plot form")
        

