# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:42:36 2024

@author: alex_liebeskind
"""
# %% =============
# Import libraries
# ================
# Loading 
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pydicom
import pypulseq as pp
from bmctool.simulate import simulate
from csaps import csaps
import argparse

# Information
# Example of calling the function with custom arguments
# eval_aptw_3t(data_flag='simulation', data_path='/path/to/data')

# data_flag:        'real_data' , 'simulation' , 're-simulation'
# data_path:        Enter path where DIDCOM data are located.
# bmsim_filename:   Enter yaml filename
# seq_filename:     enter seq filename

# Set up argparse to handle command line arguments
parser = argparse.ArgumentParser(description="EVAL_APTw_3T script")
parser.add_argument('data_flag', type=str, nargs='?', default='re_simulation',
                    help="Type of data to process: 'simulation', 're_simulation', or 'real_data'")
parser.add_argument('data_path', type=str, nargs='?', default='',
                    help="Path to the data directory")
parser.add_argument('bmsim_filename', type=str, nargs='?', default='L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml',
                    help="Filename of the BMSim configuration file")
parser.add_argument('seq_filename', type=str, nargs='?', default='APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq',
                    help="Filename of the sequence file")

args = parser.parse_args()


# Data compared
comparison_data = (
    "Glang2022/BF_volTueb1_3T_Lorentz_wo_0ppm_4cestpools_WM.yaml",
    "Heo2019_3T/WM_3T_Heo2019_4pool_bmsim.yaml",
    "Liu2013_7T/WM_3T_Liu2013_4pool_bmsim.yaml",
    "heuristic_Zaiss/WM_3T_001_bmsim.yaml",
    "vanZijl2018_3T/WM_3T_003_bmsim.yaml"
)

# Initialize lists to store data for plotting
mean_Z_corr_list = []
mean_MTRasym_list = []
offsets_list = []

# Define a colormap for distinct colors
colors = plt.cm.tab10(np.linspace(0, 1, len(comparison_data)))

# Assuming args.data_flag, args.data_path, and args.seq_filename are defined
data_flag = args.data_flag
data_path = args.data_path
seq_filename = args.seq_filename

# Define seq, config and dicom name
seq_name = Path(seq_filename)
seq_path = Path.cwd().parent / "seq-library" / seq_name.stem / seq_name
print(seq_path)
assert seq_path.is_file(), "seq file not found"

# Read in associated seq file from Pulseq-CEST library
seq = pp.Sequence()
seq.read(seq_path)

# Loop through each simulation configuration
for i, sim_location in enumerate(comparison_data):
    bmsim_filename = sim_location

    m0_offset = seq.get_definition("M0_offset")
    offsets = seq.get_definition("offsets_ppm")
    n_meas = len(offsets)

    # Re-simulate based on the configuration file
    config_name = bmsim_filename
    config_path = Path.cwd().parent / "sim-library" / config_name
    sim = simulate(config_file=config_path, seq_file=seq_path)  # Implement this function as per your library
    m_z = sim.get_zspec()[1]
    m_z = np.expand_dims(m_z, axis=1)

    # Evaluation
    M0_idx = np.where(abs(offsets) >= abs(m0_offset))[0]
    if len(M0_idx) > 0:
        M0 = np.mean(m_z[M0_idx, :], 0)
        offsets = np.delete(offsets, M0_idx)
        m_z = np.delete(m_z, M0_idx, axis=0)
        Z = m_z / M0  # Normalization
    else:
        print("m0_offset not found in offset")

    # helper function to evaluate piecewise polynomial
    def ppval(p, x):
        if callable(p):
            return p(x)
        else:
            n = len(p) - 1
            result = np.zeros_like(x)
            for i in range(n, -1, -1):
                result = result * x + p[i]
            return result
        
    # Smoothing spline interpolation
    Z_corr = np.zeros_like(Z)
    w = offsets
    dB0_stack = np.zeros(Z.shape[1])
    for ii in range(Z.shape[1]):
        if np.all(np.isfinite(Z[:, ii])):
            pp = csaps(w, Z[:, ii], smooth=0.95)
            w_fine = np.arange(-1, 1.005, 0.005)
            z_fine = ppval(pp, w_fine)

            min_idx = np.argmin(z_fine)
            dB0_stack[ii] = w_fine[min_idx]

            Z_corr[:, ii] = ppval(pp, w + dB0_stack[ii])

    # Calculate MTRasym-Spectrum
    Z_ref = Z_corr[::-1, :]
    MTRasym = Z_ref - Z_corr

    # Store mean Z_corr and MTRasym for plotting
    mean_Z_corr_list.append(np.mean(Z_corr, axis=1))
    mean_MTRasym_list.append(np.mean(MTRasym, axis=1))
    offsets_list.append(w)

# Plot Mean Z-spectrum
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
for j in range(len(comparison_data)):  # Plot all lines up to the current index
    if j == 3:  # Index 3 corresponds to "heuristic_Zaiss/WM_3T_001_bmsim.yaml"
        plt.plot(offsets_list[j], mean_Z_corr_list[j], '.-', color='black')
    else:
        plt.plot(offsets_list[j], mean_Z_corr_list[j], '.-', color=colors[j])
plt.xlim([-3.5, 3.5])
plt.gca().invert_xaxis()
plt.title("Mean Z-spectrum")
plt.xlabel("Offset (ppm)")
plt.ylabel("Normalized Intensity")

# Plot Mean MTRasym-spectrum
plt.subplot(1, 2, 2)
for j in range(len(comparison_data)):  # Plot all lines up to the current index
    if j == 3:  # Index 3 corresponds to "heuristic_Zaiss/WM_3T_001_bmsim.yaml"
        plt.plot(offsets_list[j], mean_MTRasym_list[j], '.-', color='black', label=comparison_data[j].split('/')[0])
    else:
        plt.plot(offsets_list[j], mean_MTRasym_list[j], '.-', color=colors[j], label=comparison_data[j].split('/')[0])
plt.xlim([0, 3.5])
plt.gca().invert_xaxis()
plt.title("Mean MTRasym-spectrum")
plt.xlabel("Offset (ppm)")
plt.ylabel("MTRasym")

plt.tight_layout()
plt.legend()

plt.show()
