# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 22:57:42 2025

@author: abhiv
"""

# In Class Lab 6
# Abhinav Vatsa
# Surface Brightness Profiles
# I wasn't present in class for this lab but Dr Besla provided me with a one day extension. 
# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt


# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from GalaxyMass import ComponentMass

# Create a center of mass object for M31
# I.e. an instance of the CenterOfMass class 
com_M31 = CenterOfMass('M31_000.txt', 2)  # 2 for the bulge particle
# Compute the center of mass position
com_position = com_M31.COM_P(0.1)  

print("Center of Mass Position:", com_position)

# Use the center of mass object to 
# store the x, y, z, positions and mass of the bulge particles
# be sure to correct for the COM position of M31
x = com_M31.x - com_position[0].value
y = com_M31.y - com_position[1].value
z = com_M31.z - com_position[2].value
m = com_M31.m #units of 1e10 (scientific notation)

# Determine the positions of the bulge particles in 
# cylindrical coordinates. 
cyl_r = np.sqrt(x**2 + y**2) #radial coordinate
cyl_theta = np.arctan2(y,x) # angular coordinate

def SurfaceDensity(r,m):
    """ Function that computes the surface mass density profile
    given an array of particle masses and radii 
     
    PARMETERS
    ---------
        r : array of `floats` - cyclindrical radius [kpc]
        m : array of `floats` - particle masses [1e10 Msun] 
    
    RETURNS
    -------
        r_annuli : array of `floats` -  radial bins for the 
            annuli that correspond to the surface mass density profile
    
        sigma: array of `floats` - surface mass density profile 
         [1e10 Msun/kpc^2] 
        
        
    """
    
    # Create an array of radii that captures the extent of the bulge
    # 95% of max range of bulge
    radii = np.arange(0.1, 0.95 * r.max(), 0.1)

    # create a mask to select particles within each radius
    # np.newaxis creates a virtual axis to make cyl_r_mag 2 dimensional
    # so that all radii can be compared simultaneously
    # a way of avoiding a loop - returns a boolean 
    enc_mask = r[:, np.newaxis] < radii

    # calculate mass of bulge particles within each annulus.  
    # relevant particles will be selected by enc_mask (i.e., *1)
    # outer particles will be ignored (i.e., *0)
    # axis =0 flattens to 1D
    m_enc = np.sum(m[:, np.newaxis] * enc_mask, axis=0)

    # use the difference between m_enc at adjacent radii 
    # to get mass in each annulus
    m_annuli = np.diff(m_enc) # one element less then m_enc
    
    
    # Surface mass density of stars in the annulus
    # mass in annulus / surface area of the annulus. 
    # This is in units of 1e10
    sigma = m_annuli / (np.pi * (radii[1:]**2 - radii[:-1]**2))
    # array starts at 0, but here starting at 1 and
    # subtracting radius that ends one index earlier.
    
    # Define the range of annuli
    # here we choose the geometric mean between adjacent radii
    r_annuli = np.sqrt(radii[1:] * radii[:-1]) 

    return r_annuli, sigma

# Define the surface mass density profile for the simulated bulge
# and the corresponding annuli
r_annuli, sigmaM31bulge = SurfaceDensity(cyl_r, m)

def sersicE(r, re, n, mtot):
    """ Function that computes the Sersic Profile for an Elliptical 
    System, assuming M/L ~ 1. As such, this function is also the 
    mass surface density profile. 
    
    PARMETERS
    ---------
        r: `float`
            Distance from the center of the galaxy (kpc)
        re: `float`
            The Effective radius (2D radius that contains 
            half the light) (kpc)
        n:  `float`
            sersic index
        mtot: `float`
            the total stellar mass (Msun)

    RETURNS
    -------
        I: `array of floats`
            the surface brightness/mass density
            profile for an elliptical in Lsun/kpc^2

    """
        # Set the mass-to-light ratio (M/L) to 1
    lum = mtot

    # Calculate the effective surface brightness
    Ie = lum / (7.2 * np.pi * re**2)

    # Decompose the Sersic profile
    a = (r / re) ** (1 / n)
    b = -7.67 * (a - 1)
    
    # Compute the surface brightness profile using the Sersic law
    I = Ie * np.exp(b)
    
    return I

# Create a mass profile object for M31
# using solution to Homework 5
M31mass = MassProfile("M31", 0) 

# Determine the Bulge mass profile
# use the annuli defined for the surface mass density profile
bulge_mass = M31mass.mass_enclosed_total(r_annuli).value

# Determine the total mass of the bulge
b_total = float(ComponentMass("C:/Users/abhiv/ASTR_400B/M31_000.txt", 3))*1e12
print(f"{b_total:2e}")

# Find the effective radius of the bulge, 
# Re encloses half of the total bulge mass

# Half the total bulge mass
b_half = b_total/2

# Find the indices where the bulge mass is larger than b_half
index = np.where(bulge_mass>b_half)

# take first index where Bulge Mass > b_half
# check : should match b_half
print(f"{bulge_mass[index][0]:.2e}")

# Define the Effective radius of the bulge
re_bulge = r_annuli[index][0]*3/4
print(re_bulge)

# Sersic Index = 4
SersicM31Bulge = sersicE(r_annuli, re_bulge, 4, b_total)
SersicM31Bulge

fig, ax = plt.subplots(figsize=(9,8))

#adjust tick label font size
label_size = 24
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

# Surface Density Profile
ax.loglog(r_annuli, sigmaM31bulge, lw=2, label = 'Sim Bulge', color = 'black')


# Sersic fit to the surface brightness Sersic fit
ax.loglog(r_annuli, SersicM31Bulge/1e10, linestyle = '-.', lw = 3, label = 'Sersic n=4', color = 'blue')



plt.xlabel('log r [ kpc]', fontsize=22)

# note the y axis units
plt.ylabel(r'log $\Sigma_{bulge}$ [$10^{10} M_\odot$ / kpc$^2$]', 
          fontsize=22)

plt.title('M31 Bulge', fontsize=22)

#set axis limits
plt.xlim(1,50)
plt.ylim(1e-5,0.1)

ax.legend(loc='best', fontsize=22)
fig.tight_layout()

#plt.savefig('Lab6.png')

# Data generated by the code
""""
MW COM Position (kpc): [-2.07  2.95 -1.45] kpc
MW COM Velocity (km/s): [ 0.94  6.32 -1.35] km / s
M31 COM Position (kpc): [-377.66  611.43 -284.64] kpc
M31 COM Velocity (km/s): [ 72.85 -72.14  49.  ] km / s
M33 COM Position (kpc): [-476.22  491.44 -412.4 ] kpc
M33 COM Velocity (km/s): [ 44.42 101.78 142.23] km / s
MW-M31 Separation: 769.098 kpc
MW-M31 Velocity: 117.738 km / s
M33-M31 Separation: 201.083 kpc
M33-M31 Velocity: 199.370 km / s
Galaxy: Milky Way
MW_000.txt
Halo Mass: 1.975 x 10^12 M☉
Disk Mass: 0.075 x 10^12 M☉
Bulge Mass: 0.01 x 10^12 M☉
Galaxy: Andromeda (M31)
M31_000.txt
Halo Mass: 1.921 x 10^12 M☉
Disk Mass: 0.12 x 10^12 M☉
Bulge Mass: 0.019 x 10^12 M☉
Galaxy: Triangulum (M33)
M33_000.txt
Halo Mass: 0.187 x 10^12 M☉
Disk Mass: 0.009 x 10^12 M☉
Bulge Mass: 0.0 x 10^12 M☉
Center of Mass Position: [-377.66  611.43 -284.64] kpc
1.900000e+10
1.12e+10
0.7866066361276137
"""