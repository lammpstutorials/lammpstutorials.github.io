#!/usr/bin/env python
# coding: utf-8

# In[1]:


import copy
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable


# In[2]:


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})


# In[3]:


fontsize = 15
font = {'family': 'sans', 'color':  'black', 'weight': 'normal', 'size': fontsize}
myblue = [0/ 255, 150/255, 177/ 255]
lightgray = [0.1, 0.1, 0.1]
darkgray = [0.9, 0.9, 0.9]


# In[4]:


kcal_to_joule = 4184
Na = 6.022e23
Atom = 1e-10


# #### Calculate the stress in Pa

# In[5]:


force = np.loadtxt('force.dat')
time, f_edge1, f_edge2 = force.T
f_total_kcal = f_edge2 - f_edge1 # kcal/mol/A
f_total_newton = f_total_kcal*kcal_to_joule/Na/Atom
Ly = 4e-9 # sheet size along y in meter 
dC = 3.4e-10 # effective diameter of the carbon atom
A = Ly*dC # surface area in m2
s_total = f_total_newton/A  # stress in Pascal


# #### Calculate the strain in %

# In[6]:


length = np.loadtxt('length.dat')
_, l_sheet = length.T
l_init = l_sheet[0] # Angstrom
strain = 100*(l_sheet-l_init)/l_init # %


# In[7]:


for mode, mygray in zip(['light', 'dark'], [lightgray, darkgray]):

    fig = plt.figure(figsize=(24, 4))
    ax1 = fig.add_subplot(131)
    plt.plot(strain, s_total/1e9, '-', color=myblue)
    divider = make_axes_locatable(ax1)
    ax1.set_xlabel(r"strain ($\%$)", fontdict=font)
    ax1.set_ylabel(r"stress (GPa)", fontdict=font)

    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)
    #plt.ylim(5, 65)
    ax1.minorticks_on()
    ax1.tick_params('both', length=10, width=2, which='major', direction='in')
    ax1.tick_params('both', length=6, width=1.4, which='minor', direction='in')
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_ticks_position('both')
    ax1.spines["top"].set_linewidth(2)
    ax1.spines["bottom"].set_linewidth(2)
    ax1.spines["left"].set_linewidth(2)
    ax1.spines["right"].set_linewidth(2)
    #ax1.set_xticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5])
    #ax1.set_yticks([0, 0.01, 0.02, 0.03])
    minor_locator_y = AutoMinorLocator(2)
    ax1.yaxis.set_minor_locator(minor_locator_y)
    minor_locator_x = AutoMinorLocator(2)
    ax1.xaxis.set_minor_locator(minor_locator_x)
    ax1.xaxis.label.set_color(mygray)
    ax1.yaxis.label.set_color(mygray)
    ax1.tick_params(axis='x', colors=mygray)
    ax1.tick_params(axis='y', colors=mygray)
    ax1.spines['left'].set_color(mygray)
    ax1.spines['top'].set_color(mygray)
    ax1.spines['bottom'].set_color(mygray)
    ax1.spines['right'].set_color(mygray)
    ax1.tick_params(axis='y', which='both', colors=mygray)
    ax1.tick_params(axis='x', which='both', colors=mygray)

    fig.tight_layout()
    if mode == 'light':
        plt.savefig('strain-stain-curve-light.png', bbox_inches = 'tight', pad_inches = 0.057, transparent=True)
        plt.savefig('../../../../docs/sphinx/source/tutorials/figures/graphene/strain-stain-curve-light.png', bbox_inches = 'tight', pad_inches = 0.057, transparent=True)
    else:
        plt.savefig('strain-stain-curve-dark.png', bbox_inches = 'tight', pad_inches = 0.057, transparent=True)
        plt.savefig('../../../../docs/sphinx/source/tutorials/figures/graphene/strain-stain-curve-dark.png', bbox_inches = 'tight', pad_inches = 0.057, transparent=True)        
    plt.show()

