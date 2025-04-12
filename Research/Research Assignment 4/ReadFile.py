# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 12:53:18 2025

@author: abhiv
"""
import numpy as np  
import astropy.units as u  # Importing the necessary modules

def Read(filename):
    """
    This function is designed to read a file containing Astrophysical data. 
    
    Parameters:
    filename (str): The name of the file to be read. The file is expected to have a specific format:
        - The first line contains a label and a time value in Myr (Mega-years).
        - The second line contains a label and an integer count (presumably the number of data entries).
        - The subsequent lines contain tabular data with named columns.
    
    Returns:
    (time, count, data)
        - time (Quantity): The time value converted to Astropy's unit system (Mega-years).
        - count (float): The numerical count extracted from the file (may represent the number of particles, measurements, etc.).
        - data (numpy structured array): A NumPy array with named columns containing the main dataset.
    """
    
    # Open the file 
    file = open(filename, 'r')
    
    # Read the first line 
    line1 = file.readline()  # Command to read the first line
    label, value = line1.split()  # Split the line into two parts: label and value
    time = float(value) * u.Myr  # Units of Mega-years
    
    # Read the second line 
    line2 = file.readline() # Automatically reads the next line
    label, value = line2.split()  # Split into label and numerical value
    count = float(value)  # This will give us the number of particles
    
    # Close the file after reading the necessary data
    file.close()
    
    # Load the remaining tabular data, skipping the first 3 lines
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)  
    
    return time, count, data  # Return the extracted time, count, and structured dataset


