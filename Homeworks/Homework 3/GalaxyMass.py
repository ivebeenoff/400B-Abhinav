# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 13:24:46 2025

@author: abhiv
"""
# GalaxyMass.py
import numpy as np
from ReadFile import Read  # Import the Read function

def ComponentMass(filename, particle_type):
    """
    Objective:
    Computes the total mass of a given galaxy component from the data file.

    Parameters:
    -----------
    filename : str
        The name of the file containing galaxy data.
    particle_type : int
        The type of particles to sum mass over:
        - 1: Halo
        - 2: Disk
        - 3: Bulge

    Returns:
    --------
    mass_total : float
        The total mass of the specified component in units of 10^12 solar masses,
        rounded to three decimal places.
    """
    
    # Read the file using Read function
    time, count, data = Read(filename)
    index = np.where(data['type']==particle_type)
    mass_list = 0.01*data['m'][index]
    mass_total = np.sum(mass_list)
    
    return np.round(mass_total, 3)  # Return rounded value

# List of filenames
galaxies = {
    "Milky Way": "MW_000.txt",
    "Andromeda (M31)": "M31_000.txt",
    "Triangulum (M33)": "M33_000.txt"
}

# Loop through each galaxy and calculate mass for each component
for galaxy, filename in galaxies.items():
    print(f"\nGalaxy: {galaxy}")
    print(filename)
    
    halo_mass = ComponentMass(filename, 1)
    disk_mass = ComponentMass(filename, 2)
    bulge_mass = ComponentMass(filename, 3)

    print(f"Halo Mass: {halo_mass} x 10^12 M☉")
    print(f"Disk Mass: {disk_mass} x 10^12 M☉")
    
    # Only print bulge mass if it exists (M33 might not have a bulge)
    #if bulge_mass > 0:
    print(f"Bulge Mass: {bulge_mass} x 10^12 M☉")


