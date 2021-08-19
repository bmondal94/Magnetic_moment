#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 14:52:36 2021

@author: bmondal
"""

import numpy as np


def ReadPoscar(poscar):
    """
    Generate data from the POSCAR file.

    Parameters
    ----------
    poscar : TYPE
        DESCRIPTION.

    Returns
    -------
    vector : Float 3D-array.
        3 Lattice vectors (as each row).
    poscar_pos : Float nD-array
        Cartesian atom position components of n atoms.
    ion : String
        Different ions as string as per POSCAR line-5.
    ionnumber : Integer
        Number of each ions.

    """
    filename = open(poscar, "r")
    lines = filename.readlines()
    filename.close()
    factor = float(lines[1].split()[0])
    a = np.asarray(lines[2].split()).astype(float)
    b = np.asarray(lines[3].split()).astype(float)
    c = np.asarray(lines[4].split()).astype(float)
    vector = np.array([a, b, c])*factor
    ion = np.array(lines[5].split(), dtype=str)
    ionnumber = np.array(lines[6].split(), dtype=int)
    poscar_ltype = lines[7].strip()
    poscar_pos = np.genfromtxt(
        poscar, skip_header=8, max_rows=np.sum(ionnumber), dtype=float)
    if poscar_ltype.upper().startswith('D'):
        poscar_pos = np.dot(poscar_pos, vector.T)

    return vector, poscar_pos, ion, ionnumber


def replace_multiplecations(a):
    """
    Replaces the '*' in MAGMOM from INCAR by repeatation of 2nd element.

    Parameters
    ----------
    a : Numpy array
        Array containing elemet with '*'.

    Returns
    -------
    a : Numpy array
        Array with repeated elements.

    """
    for i, s in enumerate(a):
        if '*' in s:
            element = s.split("*")
            a = np.delete(a, i)
            a = np.insert(a, i, [element[1]]*int(element[0]))
    return a


def ReadMagmom(fname, skpi_header=1):
    """
    Takes care of different version of MAGMOM=X Y Z in INCAR

    Parameters
    ----------
    fname : File name type
        The MAGMOM file name. 
    skpi_header : Integer, optional
        The 1st part of MAGMOM.dat file contains the information of OSZICAR
        MAGMOM. So you need to skip that part. How many lines are from OSZICAR.
        The default is 1.

    Returns
    -------
    magmom_all : Numpy float 2D-array
        Magnetic moments (3-components) of all atoms at different SCF steps collected from INCAR and OUTCAR.

    """
    magmom_incar = np.genfromtxt(
        fname, skip_header=skpi_header, max_rows=1, dtype=str)
    if 'MAGMOM=' in magmom_incar[0]:  # MAGMOM= 1 1 or MAGMOM=1 1
        tmp = magmom_incar[0].split("=")
        tmp = ['MAGMOM', '=', tmp[1]] if tmp[1] else ['MAGMOM', '=']
        magmom_incar = np.delete(magmom_incar, 0)
        magmom_incar = np.insert(magmom_incar, 0, tmp)
    elif magmom_incar[1] == '=':  # MAGMOM =1 1
        pass
    else:
        magmom_incar[1] = magmom_incar[1].split("=")[1]
        magmom_incar = np.insert(magmom_incar, 1, "=")
    if any("*" in s for s in magmom_incar):  # MAGMOM = 2*1
        magmom_incar = replace_multiplecations(magmom_incar)

    magmom_all = np.genfromtxt(fname, skip_header=skpi_header+1, dtype=str)
    if magmom_all.ndim > 1:
        magmom_all = np.insert(magmom_all, 0, magmom_incar, axis=0)
    else:
        magmom_all = np.vstack((magmom_incar, magmom_all))

    return magmom_all


def ReturnColor(val):
    """
    Generate color scheme in RGB mode

    Parameters
    ----------
    val : Float (< 1)
        How much red decreases in float point.

    Returns
    -------
    Tuple
        Final RGB color code float tuple.

    """
    c = [1.0, 0, 0]
    c[0] -= val
    return tuple(c)


