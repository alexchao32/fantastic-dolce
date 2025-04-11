# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 10:44:52 2025

@author: AChao
"""

# Compare MS2 results to MS2 split input files
# Ensure that I'm splitting enough to get results back from every spectrum
# Do this for all of the files

import time
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle 
import csv
import glob

time_list = []
start_time = time.time()
time_list.append(time.time())


# Read in the input file, and parse for MGF masses
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 input files/WW2DW-POCIS_AllFragmentation\split files v2')


name = []
precursor = []

with open('WW2DW-POCIS_ESI+_MSe_AllCompounds (1000).msp', "r") as file:    
    for line in file:
        if line.startswith("Name: "):
            name.append(line[6:].strip())
        elif line.startswith("PrecursorMZ: "):
            precursor.append(line[13:].strip())




# Read in the results file
os.chdir('L:/Lab/NERL_RTP_D589A_Quincy/Alex/WebApp/WebApp manuscripts/WebApp Intro paper/WW2DW analysis/MS2 WebApp results/20250318 All fragmentation, split files v2')

df1 = pd.read_csv('WW2DW_pos_1000_CFMID_results_pos.csv')
df2 = df1.drop_duplicates(subset=['ID', 'MASS_MGF'], keep='first')




















print ("\n\nTotal runtime: --- %s seconds ---" % round((time.time() - start_time), 1))  