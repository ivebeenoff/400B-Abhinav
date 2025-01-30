# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 13:26:44 2025

@author: abhiv
"""


from ReadFile import Read
import numpy as np
import astropy.units as u

def Particle_info(filename, particle_type, particle_number):
    """
    Extracts and computes fundamental properties of a given type of particle 
    from a dataset, including mass, distance, and velocity.

    Parameters:
    filename (str): Path to the data file.
    particle_type (int/str): Identifier for the type of particle (e.g., gas, stars, dark matter).
    particle_number (int): Number of the particle in the dataset.

    Returns:
    tuple: Arrays containing mass (in solar masses), distance (in kiloparsecs), 
           and velocity (in km/s) for the selected particles.
    """

    # Read the data file and extract time, total particle count, and structured data array.
    time, count, data = Read(filename)

    # Extract the type array from the dataset.
    Type = data["type"]

    # Identify indices of all particles matching the specified type.
    index = np.where(data["type"] == particle_type)

    # Convert particle masses to solar masses and apply scientific notation.
    mass_new = data["m"][index] * 10**10 * u.solMass

    # Compute the 3D spatial distance of particles from the origin (assumed center of the system).
    distance_new = np.sqrt(data["x"][index]**2 + data["y"][index]**2 + data["z"][index]**2) * u.kpc
    distance_new_ly = distance_new.to(u.lyr)  # Convert distance to light-years

    # Compute the 3D velocity magnitude of the particles.
    velocity_new = np.sqrt(data["vx"][index]**2 + data["vy"][index]**2 + data["vz"][index]**2) * u.km / u.s

    # Round computed values to three decimal places for clarity and precision
    distance_new = np.around(distance_new, 3)  # Round 3D distance (in kpc)
    distance_new_ly = np.around(distance_new_ly, 3)  # Round 3D distance (in light-years)
    velocity_new = np.around(velocity_new, 3)  # Round velocity (in km/s)
    mass_new = np.around(mass_new, 3)  # Round mass (in solar masses)

    # Extract the properties of the specified particle by index
    Mass = mass_new[particle_number-1]  # Mass of the selected particle
    Distance = distance_new[particle_number-1]  # 3D distance in kiloparsecs
    Distancely = distance_new_ly[particle_number-1]  # 3D distance in light-years
    Velocity = velocity_new[particle_number-1]  # 3D velocity in km/s

    # Return the extracted properties
    return Mass, Velocity, Distance, Distancely

# Call the function for the 100th disk particle (type 2) in the dataset
Mass, Velocity, Distance, Distancely = Particle_info("MW_000.txt", 2, 100)

# Print the values
print(Mass)  # Mass in solar masses
print(Distance)  # Distance in kiloparsecs
print(Velocity)  # Velocity in km/s
print(Distancely)  # Distance in light-years



