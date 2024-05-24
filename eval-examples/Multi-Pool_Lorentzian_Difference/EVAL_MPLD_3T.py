# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 17:42:36 2024

@author: schuerjn
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
parser.add_argument('data_flag', type=str, nargs='?', default='simulation',
                    help="Type of data to process: 'simulation', 're_simulation', or 'real_data'")
parser.add_argument('data_path', type=str, nargs='?', default='',
                    help="Path to the data directory")
parser.add_argument('bmsim_filename', type=str, nargs='?', default='WM_3T_default_7pool_bmsim.yaml',
                    help="Filename of the BMSim configuration file")
parser.add_argument('seq_filename', type=str, nargs='?', default='MultiPool_3T_002_0p9uT_80Gauss_DC50_3200ms_deepCEST.seq',
                    help="Filename of the sequence file")

args = parser.parse_args()

# Use the arguments
data_flag = args.data_flag
data_path = args.data_path
bmsim_filename = args.bmsim_filename
seq_filename = args.seq_filename


# Define seq, config and dicom name
seq_name = Path(seq_filename)
seq_path = Path.cwd().parent.parent / "seq-library" / seq_name.stem / seq_name
assert seq_path.is_file(), "seq file not found"

# 1) read in associated seq file from Pulseq-CEST library
seq = pp.Sequence()
seq.read(seq_path)
m0_offset = seq.get_definition("M0_offset")
offsets = seq.get_definition("offsets_ppm")
n_meas = len(offsets)


if data_flag == 'simulation':
    # 2a) Read in data from simulation in Pulseq folder
    seq_path_base=Path.cwd().parent.parent / "seq-library" / seq_name.stem 
    m_z = np.loadtxt(os.path.join(seq_path_base, f'M_z_{seq_name}.txt'))
    m_z = np.expand_dims(m_z, axis=1)
    
elif data_flag == 're_simulation':
    # 2b) Re-simulate
    # Implement the re-simulation using the appropriate Python library and function
    config_name = bmsim_filename
    m_z = None  # Placeholder for re-simulated data
    config_path = Path.cwd().parent.parent / "sim-library" / config_name
    sim = simulate(config_file=config_path, seq_file=seq_path)   
    m_z = sim.get_zspec()[1]
    m_z = np.expand_dims(m_z, axis=1)
elif data_flag == 'real_data':
    # 2c) Read data from measurement (DICOM)
    if data_path == '':
        dcmpath = input('Enter the path to your DICOM directory: ')
        os.chdir(dcmpath)
    else:
        dcmpath = data_path
        os.chdir(dcmpath)
    
    
    question = input('Are the DICOM Files acquired with the same protocol parameters from the PulseqCEST Library? [y/n]: ')
    if question.lower() != 'y':
        seqfile = input('Please enter the path to your seq file: ')
        seq.read(seqfile)
        offsets = seq.get_definition('offsets_ppm')
        n_meas = len(offsets)


    #read data from dicom directory
    collection = [pydicom.dcmread(os.path.join(dcmpath, filename)) for filename in sorted(os.listdir(dcmpath))]
    # extract the volume data
    V = np.stack([dcm.pixel_array for dcm in collection])
    V = np.transpose(V, (1, 2, 0))
    sz = V.shape
    V = np.reshape(V, [sz[0], sz[1], n_meas, sz[2] // n_meas]).transpose(0, 1, 3, 2)

    # Vectorization
    mask = np.squeeze(V[:, :, :, 0]) > 100
    mask_idx = np.where(mask.ravel())[0]
    V_m_z = V.reshape(-1, n_meas).T
    m_z = V_m_z[:, mask_idx]


# %% ==========
# 3) Evaluation
# =============

M0_idx = np.where(abs(offsets) >= abs(m0_offset))[0]
if len(M0_idx) > 0:
    M0 = np.mean(m_z[M0_idx, :], 0)
    offsets = np.delete(offsets, M0_idx)
    m_z = np.delete(m_z, M0_idx, axis=0)
    Z = m_z / M0  # Normalization
else:
    print("m0_offset not found in offset")


# helper function to evaluate piecewise polinomial
def ppval(p, x):
    if callable(p):
        return p(x)
    else:
        n = len(p) - 1
        result = np.zeros_like(x)
        for i in range(n, -1, -1):
            result = result * x + p[i]
        return result


# perform the smoothing spline interpolation
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


# calc of MTRasym-Spectrum
Z_ref = Z_corr[::-1, :]
MTRasym = Z_ref - Z_corr

# Vectorization Backwards
if Z.shape[1] > 1:
    V_MTRasym = np.zeros((V_m_z.shape[0], V_m_z.shape[1]), dtype=float)
    V_MTRasym[1:, mask_idx] = MTRasym
    V_MTRasym_reshaped = V_MTRasym.reshape(
        V.shape[3], V.shape[0], V.shape[1], V.shape[2]
    ).transpose(1, 2, 3, 0)

    V_Z_corr = np.zeros((V_m_z.shape[0], V_m_z.shape[1]), dtype=float)
    V_Z_corr[1:, mask_idx] = Z_corr
    V_Z_corr_reshaped = V_Z_corr.reshape(
        V.shape[3], V.shape[0], V.shape[1], V.shape[2]
    ).transpose(1, 2, 3, 0)

# %% ==========================
# 4) Plots MEAN ZSpec and MTRasym from Phantom
# =============================

plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(w, np.mean(Z_corr, axis=1), "r.-")
plt.gca().invert_xaxis()
plt.title("Mean Z-spectrum")

plt.subplot(1, 2, 2)
plt.plot(w, np.mean(MTRasym, axis=1), "b.-")
plt.xlim([0, 4])
plt.gca().invert_xaxis()
plt.title("Mean MTRasym-spectrum")
plt.show()

# %% ==================
# 5) Plot Parametric Maps from Z(3.5 ppm) and MTRasym(3.5ppm)
# =====================
if data_flag == 'real_data':
    slice_of_interest = 5  # pick slice for Evaluation
    desired_offset = 3.5
    offset_of_interest = np.where(offsets == desired_offset)[0]  # pick offset for Evaluation
    w_offset_of_interest = w[offset_of_interest]

    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.imshow(V_Z_corr_reshaped[:, :, slice_of_interest, offset_of_interest], vmin=0.5, vmax=1)
    plt.colorbar()
    plt.title("Z(Δω) = %.2f ppm" % w_offset_of_interest)
    plt.subplot(1, 2, 2)
    plt.imshow(V_MTRasym_reshaped[:, :, slice_of_interest, offset_of_interest],vmin=-0.05,vmax=0.05)
    plt.colorbar()
    plt.title("MTRasym(Δω) = %.2f ppm" % w_offset_of_interest)
    plt.show()
