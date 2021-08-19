#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 12:20:36 2021

@author: bmondal
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.animation import FuncAnimation
import subprocess
from magnetic_moment_functions import *
import argparse
import os
import sys


#%%
Q = None
def update2(ii, xs, ys, zs, u, v, w, ax, fig):
    global Q
    Q.remove()
    Q = ax.quiver(xs, ys, zs, u[ii], v[ii], w[ii], color='r', lw=2 , arrow_length_ratio=0.2)
    ax.set_title(f'Frame = {ii:03d}')
    return Q,

def onClick(event):
    global anim_running
    if anim_running:
        ani.event_source.stop()
        anim_running = False
    else:
        ani.event_source.start()
        anim_running = True
        
#%%

if __name__ == '__main__':
    
    ## ------ Setting parameters for plot -------------------------------------
    params = {'legend.fontsize': 18,
              'figure.figsize': (10,8),
             'axes.labelsize': 24,
             'axes.titlesize': 24,
             'xtick.labelsize': 24,
             'ytick.labelsize': 24,
             'errorbar.capsize':24}
    plt.rcParams.update(params)
    
    ## ------ Setting global print precision ----------------------------------
    #np.set_printoptions(precision=3, suppress=True)

    ## ------- Call parser ----------------------------------------------------
    parser = argparse.ArgumentParser(prog='magnetic_moment.py', description='This script plots the magnetic moments using VASP INCAR, OUTCAR and POSCAR file', epilog='Have fun!')
    parser.add_argument('-d', metavar='DIRECTORYNAME', default=".", help='The file path for VASP files (default: current directory). e.g. /home/mondal/VASP/test/')
    parser.add_argument('-M', action='store_true', default=False, help='Regenerate MAGMOM.dat file = True or False (default: False)')
    parser.add_argument('-Nstart', type=int, default=0, help='From which atom position you want to plot (default: 0)')
    parser.add_argument('-Nend', type=int, default=None, help='Upto which atom position you want to plot (default: till end)')
    parser.add_argument('-N', type=int, default=None, help='If instead you want to pass which group (starting with 0 in POSCAR 5th line) index of atoms for plotting (default: None)')
    parser.add_argument('-AniOSZICAR', action='store_true', default=False, help='Animate OSZICAR magnetic moments = True or False (default: False)')
    parser.add_argument('-AniFull', action='store_true', default=False, help='Animate magnetic moments from OUTCAR = True or False (default: False). If -AniFull is True then -AniOSZICAR has no effect.')
    parser.add_argument('-SaveMovie', action='store_true', default=False, help='If you want to save the Movie (default: False)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    try:
    	args = parser.parse_args()
    except:
    	#parser.print_help()
    	sys.exit(0)
    
    ## ------- Define the parser variables ------------------------------------
    dirname = args.d
    Recreate_magmom_file = args.M
    animation_oszicar = args.AniOSZICAR
    animation_full = args.AniFull
    atom_set = args.N
    start = args.Nstart
    end = args.Nend
    savefig = args.SaveMovie
    
    
    #dirname = '/home/bmondal/clusterf/project_DFT_substrate_effect/band_diagram/U3S5/KSPACING/Normal/'
    #%% -------------- CREATE MAGMOM FILE -------------------------------------
    if Recreate_magmom_file:
        print ("\n* Calling bash MagCreate.")
        subprocess.call(["/home/bmondal/script/MagCreate", dirname])
        print ("\n* MAGMOM.dat create finished.\n\n")
    
    #%% --------------- READ POSCAR FILE --------------------------------------
    vector, poscar_pos, ion, ionnumber = ReadPoscar(dirname+"POSCAR")
    
    ## ---------------- PRINT LATTICE VECTORS ---------------------------------
    with np.printoptions(precision=3, suppress=True):
        print("Lattice vectors: \n", vector)
    # print(np.array2string(vector, formatter={'float_kind':'{0:.3f}'.format}))
    
    ## --------------- Which atoms will be plotted? ---------------------------
    if atom_set:
        atom_start = ionnumber[atom_set-1] if atom_set > 0 else 0
        start = 0+atom_start; end = start+ionnumber[atom_set] # Plots 2nd set of atoms
    
    poscar_pos = poscar_pos[start:end]
    xs, ys, zs = poscar_pos[:,0], poscar_pos[:,1], poscar_pos[:,2]
    
    ## ---------------- Center of Box -----------------------------------------
    COM_box = np.mean(poscar_pos.T, axis=1)
    COM_box[2] = np.amax(vector[:,2]) 
    origin = [0,0,0]
    
    
    #%% -------------- READ MAGMOM FILE ---------------------------------------
    fname = dirname+"MAGMOM.dat"
    magmom_oszicar = np.genfromtxt(fname, comments='M', dtype=float)
    skpi_header = magmom_oszicar.size//3
    try:
        magmom_all = np.genfromtxt(dirname+"MAGMOM.dat", skip_header=skpi_header, dtype=str)
    except ValueError:
        magmom_all = ReadMagmom(fname, skpi_header)
   
    #%%
    u, v, w = {}, {}, {}
    for I in range(len(magmom_all)):
        mm = magmom_all[I][2:]
        magmom = np.split(mm, np.sum(ionnumber))
        magmom = np.array(magmom[start:end], dtype=float)
        norm_magmom = magmom #/np.amax(np.linalg.norm(magmom, axis=1))
        u[I] = norm_magmom[:,0]
        v[I] = norm_magmom[:,1]
        w[I] = norm_magmom[:,2]
 
    #%% --------------- PLOTTING ----------------------------------------------
    ## -------------- Animate MAGNETIC MOMENTS FROM OUTCAR -------------------
    if animation_full:
        fig = plt.figure()    
        ax = fig.gca(projection='3d')
        ax.scatter(xs, ys, zs, marker='o')
        ax.scatter(COM_box[0], COM_box[1], COM_box[2], marker='o', color='m', s=50)
        ax.quiver(origin, origin, origin, vector[0], vector[1], vector[2], color='k', arrow_length_ratio=0.05)
        ax.quiver(xs, ys, zs, u[0], v[0], w[0], color='m', arrow_length_ratio=0.2)
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Yaxis')
        ax.set_zlabel('Z Axis')
        anim_running = True   
        fig.canvas.mpl_connect('key_press_event', onClick)
        fig.canvas.mpl_connect('button_press_event', onClick)
        ani = FuncAnimation(fig, update2, frames=10, fargs=(xs, ys, zs, u, v, w, ax, fig) \
                            , interval=500) #, repeat=True, blit=False,repeat_delay=1))
        
    ## -------------- Animate MAGNETIC MOMENTS FROM OSZICAR -------------------
    elif animation_oszicar and skpi_header>1:
        norm_oszicar = np.linalg.norm(magmom_oszicar, axis=1)
        magmom_oszicar /= np.reshape(norm_oszicar, (len(norm_oszicar), 1))
        u, v, w = magmom_oszicar[:,0], magmom_oszicar[:,1], magmom_oszicar[:,2]
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.quiver(origin, origin, origin, [1,0,0], [0,1,0], [0,0,1], color='k', arrow_length_ratio=0.08)
        Q = ax.quiver(0, 0, 0, magmom_oszicar[0][0], magmom_oszicar[0][1], magmom_oszicar[0][2], \
                      color='r', linestyle='--', lw=3 , arrow_length_ratio=0.05)
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
        ax.set_zlim(0,1)
        plt.axis('off')
        
        anim_running = True
        fig.canvas.mpl_connect('key_press_event', onClick)
        fig.canvas.mpl_connect('button_press_event', onClick)      
        ani = FuncAnimation(fig, update2, frames=len(magmom_oszicar), fargs=(0, 0, 0, u, v, w, ax, fig), \
                            interval=500) #, repeat=True, blit=False,repeat_delay=1))
        
    else:
        ## ------------------ PLOT MAGNETIC MOMENTS from OUTCAR ---------------
        fig = plt.figure()    
        #ax = fig.add_subplot(projection='3d')
        ax = fig.gca(projection='3d')
        
        ax.scatter(xs, ys, zs, marker='o')
        ax.scatter(COM_box[0], COM_box[1], COM_box[2], marker='o', color='m', s=50)
        ax.quiver(origin, origin, origin, vector[0], vector[1], vector[2], color='k', arrow_length_ratio=0.05)
        rr = len(u)
        for I in range(rr):
            c = ReturnColor(I/rr )
            resultant_vecx, resultant_vecy, resultant_vecz = np.mean(u[I]), np.mean(v[I]), np.mean(w[I])
            print("Resultant vector-{N:d}: {X:.3f}, {Y:.3f}, {Z:.3f}".format(N=I,X=resultant_vecx, Y=resultant_vecy, \
                                                                             Z=resultant_vecz))
            ax.quiver(xs, ys, zs, u[I], v[I], w[I], color=c) 
            ax.quiver(COM_box[0], COM_box[1], COM_box[2], resultant_vecx, resultant_vecy, resultant_vecz, \
                      color=c, linestyle='--', lw=3)
        
        ax.set_xlabel('X Axis')
        ax.set_ylabel('Y Yaxis')
        ax.set_zlabel('Z Axis')
        
        ## -------------- PLOT MAGNETIC MOMENTS FROM OSZICAR ------------------
        if skpi_header>1:
            norm_oszicar = np.linalg.norm(magmom_oszicar, axis=1)
            magmom_oszicar /= np.reshape(norm_oszicar, (len(norm_oszicar), 1))
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            ax.quiver(origin, origin, origin, [1,0,0], [0,1,0], [0,0,1], color='k', arrow_length_ratio=0.08)
            for I in range(skpi_header):
                c = ReturnColor(I/rr )
                ax.quiver(0, 0, 0, magmom_oszicar[I][0], magmom_oszicar[I][1], magmom_oszicar[I][2], \
                          color=c, linestyle='--', lw=3 , arrow_length_ratio=0.05)
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            ax.set_zlim(0,1)
            plt.axis('off')
    
    ## -------------- SAVE MOVIE ----------------------------------------------
    if savefig:
        print('*** Saving movie')
        filname = dirname+'/MOVIE.mp4'
        ani.save(filname, fps=2, dpi=300)
    else:        
        plt.show()
        