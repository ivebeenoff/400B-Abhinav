#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Homework 4
# Center of Mass Position and Velocity
# Abhinav Vatsa


# In[ ]:


# I am using the template instead of writing a code from scratch as this seems more efficient. 


# In[1]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[3]:


import numpy as np
import astropy.units as u
from ReadFile import Read

class CenterOfMass:
    # Class to define COM position and velocity properties of a given galaxy and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        # create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

    def COMdefine(self, a, b, c, m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        a_com = np.sum(a * m) / np.sum(m)
        b_com = np.sum(b * m) / np.sum(m)
        c_com = np.sum(c * m) / np.sum(m)
        
        return a_com, b_com, c_com

    def COM_P(self, delta=0.1, volDec=2.0):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # First estimate
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)

        # Iterative refinement
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)

        r_max = max(r_new) / 2.0
        change = 1000.0

        while (change > delta):
            
            # Select particles within the reduced radius,
            x_new = self.x - x_COM
            y_new = self.y - y_COM
            z_new = self.z - z_COM
            r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)
            
            index2 = np.where(r_new < r_max)
            
            # Selectively retrieve those particles
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]
            
            # Calculate the COM with these particles
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2, y2, z2, m2)
            r_COM2 = np.sqrt(x_COM2**2 + y_COM2**2 + z_COM2**2)
            
            # Calculate the change in COM
            change = np.abs(r_COM - r_COM2)
            
            # Reduce r_max 
            r_max /= volDec
            
            # Reset the COM to the newly calculated values
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2
        
        # Upon reaching convergence:
        p_COM = np.array([x_COM, y_COM, z_COM]) * u.kpc
        # Rounding
        p_COM = np.round(p_COM, 2)
        
        return p_COM


    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        rv_max = 15.0 * u.kpc

        xV = self.x * u.kpc - x_COM
        yV = self.y * u.kpc - y_COM
        zV = self.z * u.kpc - z_COM
        rV = np.sqrt(xV**2 + yV**2 + zV**2)

        indexV = np.where(rV < rv_max)

        vx_new = self.vx[indexV]
        vy_new = self.vy[indexV]
        vz_new = self.vz[indexV]
        m_new = self.m[indexV]

        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)

        v_COM = np.array([vx_COM, vy_COM, vz_COM]) * u.km / u.s
        return np.round(v_COM, 2)


# In[5]:


# Create a Center of Mass object for the Milky Way (MW)
MW_COM = CenterOfMass("MW_000.txt", 2)

# Create a Center of Mass object for Andromeda (M31)
M31_COM = CenterOfMass("M31_000.txt", 2)

# Create a Center of Mass object for Triangulum (M33)
M33_COM = CenterOfMass("M33_000.txt", 2)


# In[6]:


# 1.
# MW: Store the position and velocity COM
MW_COM_p = MW_COM.COM_P(0.1)
print("MW COM Position (kpc):", MW_COM_p)

MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
print("MW COM Velocity (km/s):", MW_COM_v)

# M31: Store the position and velocity COM
M31_COM_p = M31_COM.COM_P(0.1)
print("M31 COM Position (kpc):", M31_COM_p)

M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
print("M31 COM Velocity (km/s):", M31_COM_v)

# M33: Store the position and velocity COM
M33_COM_p = M33_COM.COM_P(0.1)
print("M33 COM Position (kpc):", M33_COM_p)

M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
print("M33 COM Velocity (km/s):", M33_COM_v)


# In[11]:


# 2. 
# Separation and Velocity Between MW and M31
# Compute separation between MW and M31
MW_M31_sep = np.sqrt((MW_COM_p[0] - M31_COM_p[0])**2 +
                      (MW_COM_p[1] - M31_COM_p[1])**2 +
                      (MW_COM_p[2] - M31_COM_p[2])**2)

# Compute relative velocity
MW_M31_vel = np.sqrt((MW_COM_v[0] - M31_COM_v[0])**2 +
                      (MW_COM_v[1] - M31_COM_v[1])**2 +
                      (MW_COM_v[2] - M31_COM_v[2])**2)

# Print results rounded to three decimal places
print(f"MW-M31 Separation: {MW_M31_sep:.3f}")
print(f"MW-M31 Velocity: {MW_M31_vel:.3f}")

# 3. 
# Compute separation between M33 and M31
M33_M31_sep = np.sqrt((M33_COM_p[0] - M31_COM_p[0])**2 +
                      (M33_COM_p[1] - M31_COM_p[1])**2 +
                      (M33_COM_p[2] - M31_COM_p[2])**2)

# Compute relative velocity
M33_M31_vel = np.sqrt((M33_COM_v[0] - M31_COM_v[0])**2 +
                      (M33_COM_v[1] - M31_COM_v[1])**2 +
                      (M33_COM_v[2] - M31_COM_v[2])**2)

# Print results rounded to three decimal places
print(f"M33-M31 Separation: {M33_M31_sep:.3f}")
print(f"M33-M31 Velocity: {M33_M31_vel:.3f}")





# In[ ]:


# 4
# The iterative process for finding the COM is essential because MW and M31 are in the process of merging, which means their mass distribution isn’t uniform or stable. 
# There are tidal tails, streams, and other particles that aren’t tightly bound, and if we just took a simple average, those could throw off our results. 
# By shrinking the search radius step by step, we can focus on the densest, most gravitationally bound part of the galaxy and get a more accurate COM.
# This methodology is useful when dealing with galaxies that are interacting and not in perfect equilibrium, hence it was invoked in this instance.

