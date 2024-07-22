# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 16:27:34 2022

@author: user
"""

from mesh_in_MPM_2D import BC
import numpy as np
import matplotlib.pyplot as plt
from sys import exit

FilePath = r"particle_output.txt"

ptc_nums = [16]
node_num = 169

StartTime = 0.0
EndTime = 0.5
dt = 1e-4
MPMRecStep = 1

N = int(round((EndTime/dt)/MPMRecStep*ptc_nums[0]+ptc_nums[0], 8))
xpA = [0]*N
with open(FilePath, 'r') as file:
    i, r = 0, 0
    for text in file:
        if (r != 0):
            xpA[i] = text.split()
            xpA[i][0] = float(xpA[i][2])
            xpA[i][1] = float(xpA[i][3])
            xpA[i] = np.array([xpA[i][0] ,xpA[i][1]])
            i += 1
        r += 1

xn = [0]*node_num    
with open(r"nodal_position.txt", 'r') as file:
    i = 0
    for text in file:
        xn[i] = text.split()
        xn[i][0] = float(xn[i][0])
        xn[i][1] = float(xn[i][1])
        xn[i] = np.array(xn[i])
        i += 1

Min_bc_val = np.min(xn)
Max_bc_val = np.max(xn)

bc1 = BC("slip", p1=np.array([-1.5, -1.5]), p2=np.array([1.5, -1.5]))
bc2 = BC("slip", p1=np.array([1.5,  -1.5]), p2=np.array([1.5,  1.5]))
bc3 = BC("slip", p1=np.array([1.5,   1.5]), p2=np.array([-1.5, 1.5]))
bc4 = BC("slip", p1=np.array([-1.5,  1.5]), p2=np.array([-1.5,-1.5]))
bc1.set_normal(nbc=np.array([0.0,  1.0]))
bc2.set_normal(nbc=np.array([-1.0, 0.0]))
bc3.set_normal(nbc=np.array([0.0, -1.0]))
bc4.set_normal(nbc=np.array([1.0,  0.0]))
bc_array = [bc1, bc2, bc3, bc4]
       
def plot_fig(i, k):
    T90 = np.array([[0, -1],
                    [1, 0]])

    fig,ax = plt.subplots(figsize=(16,9), dpi=150)
    plt.ion()
    for j in range(node_num-1):
        if (xn[j][1] == xn[j+1][1]):
            xno = np.array([xn[j][0], xn[j][1]]) 
            xnl = np.array([xn[j+1][0], xn[j+1][1]])
            ax.plot([xno[0], xnl[0]], [xno[1], xnl[1]], "-", lw=0.3, color="tab:orange")
            
            xn2 = T90.dot(xno)
            xn3 = T90.dot(xnl)
            ax.plot([xn2[0], xn3[0]], [xn2[1], xn3[1]], "-", lw=0.3, color="tab:orange")

    for j in range(int(i*ptc_nums[0]), int((i+1)*ptc_nums[0])):
        ax.plot(xpA[j][0], xpA[j][1], "o", markersize=3, color="tab:blue")
    
    ax.plot(1e+2, 1e+2, "o", markersize=10, color="tab:blue", label="block") 
    ax.plot(1e+3, 1e+3, "s", markersize=10, color="tab:orange", label="Grid", markerfacecolor='none')
    ax.plot(1e+2, 1e+2, "-", lw=7, color="gray", label="BC")
    
    for bc in bc_array:
        poo = bc.p1
        pll = bc.p2
        ax.plot([poo[0], pll[0]], [poo[1], pll[1]], "-", lw=5, color="gray")
    
    ax.set_xlim(Min_bc_val-0.2, Max_bc_val+0.2)
    ax.set_ylim(Min_bc_val-0.2, Max_bc_val+0.2)
    
    ax.set_xticks(np.linspace(Min_bc_val, Max_bc_val, 5))
    ax.set_yticks(np.linspace(Min_bc_val, Max_bc_val, 5))
    
    ax.set_xlabel("x(m)", fontsize=25)
    ax.set_ylabel("y(m)", fontsize=25)

    ax.set_title("t = %.2f (s)"%tns[i], fontsize=25)
    ax.set_aspect('equal', 'box')
    ax.tick_params(labelsize=25)
    ax.ticklabel_format(style='sci', scilimits=(-0,0), axis='y')
    # ax.legend(fontsize=25, loc=9, framealpha=1, ncol=2)
    ax.legend(bbox_to_anchor=(1.2, 1), loc=9, borderaxespad=0, fontsize=25, ncol=1)
    plt.ioff()
    plt.savefig(r"C:\Users\lintim0622\source\repos\MPM\temp fig\%08d_all_responses.png"%(i+1+k))
    plt.close("all")

plt.rcParams["font.family"] = "Times New Roman"
tns = np.arange(StartTime, EndTime+dt*MPMRecStep/10.0, dt*MPMRecStep)
for i in range(int(N/ptc_nums[0])):
    if (i % 100 == 0):
        plot_fig(i, k=0)
        
    if (i == (int(N/ptc_nums[0]))):
        for k in range(5):
            plot_fig(i, k)
