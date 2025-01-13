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
# data_path:        Enter path where DICOM data are located.
# bmsim_filename:   Enter yaml filename
# seq_filename:     enter seq filename

# Set up argparse to handle command line arguments
parser = argparse.ArgumentParser(description="EVAL_APTw_3T script")
parser.add_argument('data_flag', type=str, nargs='?', default='real_data',
                    help="Type of data to process: 'simulation', 're_simulation', or 'real_data'")
parser.add_argument('data_path', type=str, nargs='?', default='',
                    help="Path to the data directory")
parser.add_argument('bmsim_filename', type=str, nargs='?', default='/Users/alexliebeskind/Desktop/pulseq-cest-library/sim-library/phantoms/l-arginin/WM_3T_default_2pool_bmsim_jschuere_water_cest_zyste_20mM.yaml',
                    help="Filename of the BMSim configuration file")
parser.add_argument('seq_filename', type=str, nargs='?', default='APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq',
                    help="Filename of the sequence file")
parser.add_argument('interpolation', type=str, nargs='?', default='spline',
                    help="Type of interpolation (linear or spline)")

args = parser.parse_args()

# Use the arguments
data_flag = args.data_flag
data_path = args.data_path
bmsim_filename = args.bmsim_filename
seq_filename = args.seq_filename
interpolation = args.interpolation

# Define seq, config and dicom name
seq_name = Path(seq_filename)
seq_path = Path.cwd() / "seq-library" / seq_name.stem / seq_name
assert seq_path.is_file(), "seq file not found"
ROI = 'n'
x_min = 0
x_max = 0
y_min = 0
y_max = 0

# 1) read in associated seq file from Pulseq-CEST library
seq = pp.Sequence()
seq.read(seq_path)
m0_offset = seq.get_definition("M0_offset")
offsets = seq.get_definition("offsets_ppm")
n_meas = len(offsets)

if data_flag == 'simulation':
    # 2a) Read in data from simulation in Pulseq folder
    seq_path_base=Path.cwd().parent.parent / "Desktop" / "pulseq-cest-library" / "seq-library" / seq_name.stem 
    m_z = np.loadtxt(os.path.join(seq_path_base, f'M_z_{seq_name}.txt'))
    m_z = np.expand_dims(m_z, axis=1)

    # Plot m_z
    plt.figure(figsize=(5, 4))
    plt.plot(np.mean(m_z, axis=1), ".-")
    plt.xlabel(r'$\Delta\omega$ [offset]')
    plt.ylabel(r'Z($\Delta\omega$)')
    plt.title("Z-spectrum")
    plt.show()
    
elif data_flag == 're_simulation':
    # 2b) Re-simulate
    # Implement the re-simulation using the appropriate Python library and function
    config_name = bmsim_filename
    m_z = None  # Placeholder for re-simulated data
    config_path = Path.cwd().parent.parent / "Desktop" / "pulseq-cest-library" / "sim-library" / "phantoms" / "l-arginin" / "WM_3T_default_2pool_bmsim_jschuere_water_cest_zyste_20mM.yaml"
    sim = simulate(config_file=config_path, seq_file=seq_path)   
    m_z = sim.get_zspec()[1]
    m_z = np.expand_dims(m_z, axis=1)

    # Plot m_z
    plt.figure(figsize=(5, 4))
    plt.plot(offsets, np.mean(m_z, axis=1), ".-")
    plt.xlabel(r'$\Delta\omega$ [ppm]')
    plt.ylabel(r'Z($\Delta\omega$)')
    plt.gca().invert_xaxis()
    plt.xlim([3.5, -3.5])
    plt.ylim([0, 1])
    plt.title("Z-spectrum")
    plt.show()

elif data_flag == 'real_data':
    # 2c) Read data from measurement (DICOM)
    if data_path == '':
        dcmpath = '/Users/alexliebeskind/Desktop/pulseq-cest-library/data/PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015'
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
    
    # Plot the non-normalized Z-spectrum
    # plt.figure(figsize=(8, 6))
    # plt.plot(offsets, np.mean(m_z, axis=1), "g.-", label="Non-Normalized Z-Spectrum")
    # plt.xlim([3.5, -3.5])  # Set x-axis range
    # plt.xlabel(r'$\Delta\omega$ [ppm]')
    # plt.ylabel(r'Z($\Delta\omega$) (Raw Signal)')
    # plt.title("Non-Normalized Z-Spectrum")
    # plt.grid()
    # plt.legend()
    # plt.show()

    ROI = input('Do you want to specify an ROI [y/n]: ')
    if ROI == 'y':
        slice_of_interest = 5
        x_min = 145
        x_max = 150
        y_min = 95
        y_max = 100

        # Create mask for the specified ROI
        mask_ROI = np.zeros_like(mask)
        mask_ROI[x_min:x_max, y_min:y_max, slice_of_interest] = True
        mask_idx_ROI = np.where(mask_ROI.ravel())[0]
        V_m_z_ROI = V.reshape(-1, n_meas).T
        m_z_ROI = V_m_z_ROI[:, mask_idx_ROI]

        # Plot m_z
        plt.figure(figsize=(5, 4))
        plt.plot(offsets[:len(offsets)], np.mean(m_z_ROI[:len(offsets)], axis=1), ".-")
        plt.xlim([-3.5,3.5])
        plt.gca().invert_xaxis()
        plt.xlabel(r'$\Delta\omega$ [offset]')
        plt.ylabel(r'Z($\Delta\omega$)')
        plt.title("Z-spectrum")
        plt.show()

        # Plot m_z
        # plt.figure(figsize=(5, 4))
        # plt.plot(offsets[:len(offsets)], np.mean(m_z[:len(offsets)], axis=1), ".-")
        # #plt.xlim([-3.5,3.5])
        # #plt.ylim(0,1)
        # plt.gca().invert_xaxis()
        # plt.xlabel(r'$\Delta\omega$ [offset]')
        # plt.ylabel(r'Z($\Delta\omega$)')
        # plt.title("Z-spectrum")
        # plt.show()

# Plot m_z
    # plt.figure(figsize=(5, 4))
    # plt.plot(offsets, np.mean(m_z, axis=1), ".-")
    # plt.gca().invert_xaxis()
    # plt.xlim([3.5, -3.5])
    # plt.xlabel(r'$\Delta\omega$ [offset]')
    # plt.ylabel(r'Z($\Delta\omega$)')
    # plt.title("Raw Z-spectrum")
    # plt.show()


# %% ==========
# 3) Evaluation
# =============

# Extract the ROI for slice of interest in z dimension
M0_idx = np.where(abs(offsets) >= abs(m0_offset))[0]
if len(M0_idx) >= 0:
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

from scipy.interpolate import interp1d

if interpolation == "linear":
    # perform linear interpolation
    Z_corr = np.zeros_like(Z)
    w = offsets
    dB0_stack = np.zeros(Z.shape[1])

    for ii in range(Z.shape[1]):
        if np.all(np.isfinite(Z[:, ii])):
            # Create linear interpolation function
            f = interp1d(w, Z[:, ii], kind='linear', fill_value='extrapolate')
            
            # Interpolate values at fine grid points
            w_fine = np.arange(-1, 1.0001, 0.0001)
            z_fine = f(w_fine)

            # Find index of minimum value
            min_idx = np.argmin(z_fine)
            dB0_stack[ii] = w_fine[min_idx]
            
            # modified based on 9/24 -- need to account for second peak
            # dB0_stack = dB0_stack - 0.028

            # Interpolate corrected values
            Z_corr[:, ii] = f(w + dB0_stack[ii]-0.028)

elif interpolation == "spline": 
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

            # modified based on 9/24 -- need to account for second peak
            # dB0_stack = dB0_stack - 0.028

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
# 4) Plot MEAN ZSpec and MTRasym from Phantom
# =============================

if ROI == 'n':
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(w, np.mean(Z_corr, axis=1), ".-")
    plt.gca().invert_xaxis()
    plt.xlim([3.5, -3.5])
    plt.ylim([0, 1])
    plt.xlabel(r'$\Delta\omega$ [ppm]')
    plt.ylabel(r'Z($\Delta\omega$)')
    plt.title("Z-spectrum")
    plt.subplot(1, 2, 2)
    plt.plot(w, np.mean(MTRasym, axis=1), ".-")
    plt.xlim([0, 3.5])
    plt.ylim([-0.18, 0.1])
    plt.gca().invert_xaxis()
    plt.xlabel(r'$\Delta\omega$ [ppm]')
    plt.ylabel(r'$MTR_{asym}(\Delta\omega)$')
    plt.title(r'$MTR_{asym}$')
    plt.tight_layout()
    plt.show()
else:
    # Extract the ROI for slice of interest in z dimension
    V_Z_data = V_Z_corr_reshaped[x_min:x_max, y_min:y_max, slice_of_interest, :]
    V_MTRasym_data = V_MTRasym_reshaped[x_min:x_max, y_min:y_max, slice_of_interest, :]

    # Calculate the average spectrum across the ROI for each w
    Z_spectrum = np.mean(V_Z_data, axis=(0, 1))
    V_MTRasym_spectrum = np.mean(V_MTRasym_data, axis=(0, 1))

    Z_spectrum = Z_spectrum[:len(w)]
    V_MTRasym_spectrum = V_MTRasym_spectrum[:len(w)]


    #------------
    if interpolation == "linear":
        # Perform linear interpolation
        f = interp1d(w, Z_spectrum, kind='linear', fill_value='extrapolate')
        w_fine = np.arange(-1, 1.0005, 0.0005)
        z_fine = f(w_fine)
        min_idx = np.argmin(z_fine)
        dB0 = w_fine[min_idx]
        Z_spectrum_corr = f(w + dB0)

        # Correct MTR asymmetry spectrum using the same dB0 shift
        f_mtr = interp1d(w, V_MTRasym_spectrum, kind='linear', fill_value='extrapolate')
        V_MTRasym_spectrum_corr = f_mtr(w + dB0)

    elif interpolation == "spline":
        # Perform the smoothing spline interpolation
        pp = csaps(w, Z_spectrum, smooth=0.99)
        w_fine = np.arange(-1, 1.005, 0.005)
        z_fine = ppval(pp, w_fine)
        min_idx = np.argmin(z_fine)
        dB0 = w_fine[min_idx]

        #0.028 adjustment here
        Z_spectrum_corr = ppval(pp, w + dB0)

        # Correct MTR asymmetry spectrum using the same dB0 shift
        pp_mtr = csaps(w, V_MTRasym_spectrum, smooth=0.99)
        V_MTRasym_spectrum_corr = ppval(pp_mtr, w + dB0)
    
    # Plot the average spectra
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(w, Z_spectrum_corr, ".-")
    plt.xlim([-3.5, 3.5])
    plt.ylim([0, 1])
    plt.gca().invert_xaxis()
    plt.xlabel(r'$\Delta\omega$ [ppm]')
    plt.ylabel(r'Z($\Delta\omega$)')
    plt.title("Z-spectrum")
    plt.subplot(1, 2, 2)
    plt.plot(w, V_MTRasym_spectrum_corr, ".-")
    plt.xlim([0, 3.5])
    plt.ylim([-0.18, 0.1])
    plt.gca().invert_xaxis()
    plt.xlabel(r'$\Delta\omega$ [ppm]')
    plt.ylabel(r'$MTR_{asym}(\Delta\omega)$')
    plt.title(r'$MTR_{asym}$')
    plt.tight_layout()

    plt.show()

# %% ==================
# 5) Plot Parametric Maps from Z(3.5 ppm) and MTRasym(3.5ppm)
# =====================
if data_flag == 'real_data':
    slice_of_interest = 5  # pick slice for Evaluation
    desired_offset = 3
    offset_of_interest = np.where(offsets == desired_offset)[0]  # pick offset for Evaluation
    w_offset_of_interest = w[offset_of_interest]

    plt.figure(figsize=(10, 4))
    ax1 = plt.subplot(1, 2, 1)
    plt.imshow(V_Z_corr_reshaped[:, :, slice_of_interest, offset_of_interest], vmin=0.5, vmax=1)
    plt.colorbar()
    plt.title(r'Z($\Delta\omega$) = %.2f ppm' % w_offset_of_interest)
    
    x_max, y_max = y_max, x_max
    x_min, y_min = y_min, x_min

    if ROI == 'y': 
        rect1 = plt.Rectangle((x_min, y_min), (x_max - x_min), (y_max - y_min), linewidth=1, edgecolor='red', facecolor='none')
        ax1.add_patch(rect1)

    ax2 = plt.subplot(1, 2, 2)
    plt.imshow(V_MTRasym_reshaped[:, :, slice_of_interest, offset_of_interest],vmin=-0.05,vmax=0.05)
    plt.colorbar()
    plt.title(r'$MTR_{asym}(\Delta\omega)$ = %.2f ppm' % w_offset_of_interest)

    if ROI == 'y':
        rect2 = plt.Rectangle((x_min, y_min), (x_max - x_min), (y_max - y_min), linewidth=1, edgecolor='red', facecolor='none')
        ax2.add_patch(rect2)

    plt.show()