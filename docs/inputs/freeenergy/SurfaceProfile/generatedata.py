#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import random
import copy


# ## Homemade functions
# 
# source : https://github.com/simongravelle/python-for-lammps
from numpy.linalg import norm
def neighborsearch(neighbor,molecule,cptatm, x, y, z, Lx, Ly, Lz):
    '''Search neighbor in a box and return the closest distance found'''
    box = np.array([Lx, Ly, Lz])
    minr = 10
    for m in molecule.T:
        x0 = m[0] + x
        y0 = m[1] + y
        z0 = m[2] + z
        dxdydz = np.remainder(neighbor[:cptatm].T - np.array([x0,y0,z0]) + box/2., box) - box/2.
        minr = np.min([minr,np.min(norm(dxdydz,axis=1))])
    return minr

from scipy.spatial.transform import Rotation as R
def randomorientation(XYZ):
    '''3D aleatory rotation of molecule/particule coordinate'''
    rotation_degrees = random.randint(0,9000)/100
    rotation_radians = np.radians(rotation_degrees)
    rotation_axis = np.array([random.randint(0,100)/100, random.randint(0,100)/100, random.randint(0,100)/100])
    rotation_axis /= np.linalg.norm(rotation_axis)
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    mol_rotated = rotation.apply(XYZ)    
    return mol_rotated.T

def randomlocation(Lx,Ly,Lz):
    '''Choose a random location within a given box'''
    txlo, txhi = -Lx/2, Lx/2
    tylo, tyhi = -Ly/2, Ly/2
    tzlo, tzhi = -Lz/2, Lz/2    
    x = random.randint(1,1000)/1000*(txhi-txlo)+txlo
    y = random.randint(1,1000)/1000*(tyhi-tylo)+tylo
    z = random.randint(1,1000)/1000*(tzhi-tzlo)+tzlo
    return x, y, z


# ## parameter choice
nbethanol = 10 # desired number of ethanol molecules
h = 40 # distance between the walls
layer = 12 # desired layer size


# ## NaCl layer dimension
dnacl = 2.84 # Na-Cl typical distance
nx, ny, nz = 10, 10, 4
Lx, Ly, Lz = nx*dnacl, ny*dnacl, nz*dnacl
print('The desired Na-Cl layer dimensions are '+str(Lx)+' A x ' +str(Ly)+' A x '+str(Lz)+' A')


# ## Lammps box size
txlo, txhi = -Lx/2, Lx/2
tylo, tyhi = -Ly/2, Ly/2
tzlo, tzhi = -Lz/2-h/2, Lz/2+h/2


# ## create the NaCl wall
basestructure = np.loadtxt('NaCl/Positions.dat') # import the 8 basic atoms to replicates
# replicate the initial structure
naclwall = copy.deepcopy(basestructure)
for xx in np.arange(txlo+dnacl/2,txhi,2*dnacl):
    for yy in np.arange(tylo+dnacl/2,tyhi,2*dnacl):
        for zz in np.arange(tzlo+dnacl/2,tzlo+layer,2*dnacl):
            naclwall = np.append(naclwall,basestructure+[0,0,0,0,xx,yy,zz], axis=0)
naclwall = naclwall[8:]
# renumber atoms ids
for n in range(len(naclwall)):
    naclwall[n,0] = np.int64(n+1)
# measure the effective width of the NaCl block
efflayer = np.max(naclwall.T[6]) - np.min(naclwall.T[6])
print('The effective layer width is '+str(efflayer)+' A')
# recenter the wall
naclwall.T[6] -= np.mean(naclwall.T[6])-Lz/2-h/2 


# ## create data lammps file
cptatm = 0
cptbnd = 0
cptang = 0
cptdih = 0
cptmol = 0
atoms = np.zeros((1000000,7))
bonds = np.zeros((1000000,4))
angles = np.zeros((1000000,5))
dihedrals = np.zeros((1000000,6))


# ## place the wall
cptmol += 1
for m in naclwall:
    atoms[cptatm] = m[0], cptmol, m[2], m[3], m[4], m[5], m[6]
    cptatm += 1


# ## add ethanol molecules

posEth = np.loadtxt('EthanolMolecule/Positions.dat')
posEth.T[2] += 2 # shift to account for NaCl atoms
bndEth = np.loadtxt('EthanolMolecule/Bonds.dat')
angEth = np.loadtxt('EthanolMolecule/Angles.dat')
dihEth = np.loadtxt('EthanolMolecule/Dihedrals.dat')

while cptmol-1 < nbethanol:
    x, y, z = randomlocation(Lx,Ly,Lz+h)
    posEth.T[4:7] = randomorientation(posEth.T[4:7].T)
    mind = neighborsearch(atoms.T[4:7],posEth.T[4:7],cptatm, x, y, z, Lx, Ly, Lz+h)
    if mind>3:
        cptmol += 1
        for m in bndEth:
            bonds[cptbnd] = cptbnd+1, m[1], m[2]+cptatm, m[3]+cptatm
            cptbnd += 1
        for m in angEth:
            angles[cptang] = cptang+1, m[1], m[2]+cptatm, m[3]+cptatm, m[4]+cptatm
            cptang += 1
        for m in dihEth:
            dihedrals[cptdih] = cptdih+1, m[1], m[2]+cptatm, m[3]+cptatm, m[4]+cptatm, m[5]+cptatm
            cptdih += 1
        for m in posEth:
            atoms[cptatm] = cptatm+1, cptmol, m[2], m[3], m[4]+x, m[5]+y, m[6]+z
            cptatm += 1


# ## remove excess lines

atoms = atoms[0:cptatm]       
bonds = bonds[0:cptbnd]    
angles = angles[0:cptang]
dihedrals = dihedrals[0:cptdih]


# ## write LAMMPS data file

f = open("data.lammps", "w")
f.write('# LAMMPS data file \n\n')
f.write(str(cptatm)+' atoms\n')
f.write(str(cptbnd)+' bonds\n')
f.write(str(cptang)+' angles\n')
f.write(str(cptdih)+' dihedrals\n')
f.write('\n')
f.write(str(int(7))+' atom types\n')
f.write(str(int(5))+' bond types\n')
f.write(str(int(6))+' angle types\n')
f.write(str(int(3))+' dihedral types\n')
f.write('\n')
f.write(str(txlo)+' '+str(txhi)+' xlo xhi\n')
f.write(str(tylo)+' '+str(tyhi)+' ylo yhi\n')
f.write(str(tzlo)+' '+str(tzhi)+' zlo zhi\n')
f.write('\n')
f.write('Atoms\n')
f.write('\n')
for nlin in range(len(atoms)):
    newline = atoms[nlin]
    for col in range(len(newline)):
        if col < 3:
            f.write(str(int(newline[col]))+' ')
        else :
            f.write(str(newline[col])+' ')
    f.write('\n')
f.write('\n')
f.write('Bonds\n')
f.write('\n')
for nlin in range(len(bonds)):
    newline = bonds[nlin]
    for col in range(len(newline)):
        f.write(str(int(newline[col]))+' ')
    f.write('\n')
f.write('\n')
f.write('Angles\n')
f.write('\n')
for nlin in range(len(angles)):
    newline = angles[nlin]
    for col in range(len(newline)):
        f.write(str(int(newline[col]))+' ')
    f.write('\n')
f.write('\n')
f.write('Dihedrals\n')
f.write('\n')
for nlin in range(len(dihedrals)):
    newline = dihedrals[nlin]
    for col in range(len(newline)):
        f.write(str(int(newline[col]))+' ')
    f.write('\n')
f.close()
