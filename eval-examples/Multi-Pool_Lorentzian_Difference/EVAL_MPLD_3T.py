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
parser.add_argument('data_flag', type=str, nargs='?', default='real_data',
                    help="Type of data to process: 'simulation', 're_simulation', or 'real_data'")
parser.add_argument('data_path', type=str, nargs='?', default=r'\\141.67.249.47\MRTransfer\pulseq_zero\sequences\seq240524\CEST\MultiPool_3T_002_0p9uT_80Gauss_DC50_3200ms_deepCEST_cFA_img.npy',
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
    
    path = Path(data_path)
    if path.is_file() and path.suffix == '.npy':  # read npy
        V=np.abs(np.load(data_path))
        sz = V.shape
        print(sz)
        # Vectorization
        mask = np.squeeze(V[:, :, :, 0]) > np.mean(np.abs(V))/5
        mask_idx = np.where(mask.ravel())[0]
        V_m_z = V.reshape(-1, n_meas).T
        m_z = V_m_z[:, mask_idx]
    
    else:
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
        pp = csaps(w, Z[:, ii], smooth=0.999)
        w_fine = np.arange(-1, 1.005, 0.005)
        z_fine = ppval(pp, w_fine)

        min_idx = np.argmin(z_fine)
        dB0_stack[ii] = w_fine[min_idx]

        Z_corr[:, ii] = ppval(pp, w + dB0_stack[ii])


from sklearn.decomposition import PCA

def pca_denoise(data, n_components):
    # Perform PCA
    pca = PCA(n_components=n_components)
    transformed_data = pca.fit_transform(data.T)
    
    # Inverse transform to reconstruct the data
    reconstructed_data = pca.inverse_transform(transformed_data)
    
    # Reshape back to original 4D shape
    denoised_data = reconstructed_data.T
    
    return denoised_data


print("Original Data Shape:", Z_corr.shape)

Z_corr = pca_denoise(Z_corr, 10)

print("Denoised Data Shape:", Z_corr.shape)



#%% compute the analytic function of wasabi curve    
from scipy.optimize import curve_fit

def lorentzfit4pool_rel_fixPeaks(w,*p):
    # Constants for peak positions
    dwMT = -2.5
    dwAPT = 3.5
    dwNOE = -3.5

    # Calculating the model
    z = (p[0] 
         - p[1] * p[2]**2 / 4 / (p[2]**2 / 4 + (w - p[3])**2)
         - p[4] * p[5]**2 / 4 / (p[5]**2 / 4 + (w - dwAPT - p[3])**2)
         - p[6] * p[7]**2 / 4 / (p[7]**2 / 4 + (w - dwNOE - p[3])**2)
         - p[8] * p[9]**2 / 4 / (p[9]**2 / 4 + (w - dwMT - p[3])**2))
    return z

# Initial parameters
p0 = [1, 0.90, 2.3, 0, 0.0, 2, 0.1, 4, 0.1, 60]
bounds = ([0.5, 0.40, 1, -1, 0.0, 1, 0.0, 2, 0.0025, 30], 
          [1, 1, 6, 1, 0.2, 10, 0.2, 12.5, 0.3, 100])   
 
# perform the fit 
amide_stack = np.zeros(Z.shape[1]) # we initialize an array dB0_stack of zeros with a length equal to the number of columns in the Z array
rNOE_stack = np.zeros(Z.shape[1]) # same as the previous line
MT_stack = np.zeros(Z.shape[1]) # same as the previous line
FIT_stack = np.zeros((Z.shape[1],np.size(w))) # same as the previous line
p_stack = np.zeros((Z.shape[1],np.size(p0))) # same as the previous line


for ii in range(Z.shape[1]):
    if np.all(np.isfinite(Z[:, ii])): # we check if all values in the current column ii of Z are finite (i.e., not NaN or infinity)
       
        try:
            p = curve_fit(lorentzfit4pool_rel_fixPeaks, w, Z_corr[:,ii],  p0=p0,bounds=bounds,xtol=1e-4, ftol=1e-3,maxfev=300,full_output=True, method='dogbox') # curve_fit 
            if (ii%1000): print(f'pixel {ii} of {Z.shape[1]} : {p[0]}.2f')
            amide_stack[ii]= p[0][4] # amide
            rNOE_stack[ii]= p[0][6] # amide
            MT_stack[ii]=p[0][8] # amide
            FIT_stack[ii,:]=lorentzfit4pool_rel_fixPeaks(w,*p[0])
            p_stack[ii,:]=p[0]
             
        except Exception:
           print(f'pixel {ii} of {Z.shape[1]} failed')
           amide_stack[ii]= np.nan  # in case an exception is raised during the fitting process, we set that particular column to NaN (Not-a-Number)
           rNOE_stack[ii]= np.nan
           MT_stack[ii]= np.nan
           FIT_stack[ii,:]= np.nan
           p_stack[ii,:]

if Z.shape[1] == 1:
    fig, ax = plt.subplots()
    
    fit0=lorentzfit4pool_rel_fixPeaks(w,*p0)
    fit=lorentzfit4pool_rel_fixPeaks(w,*p[0])
    data=Z_corr[:,ii]
    
    ax.plot(w, fit0, label='fit')
    ax.plot(w, fit, label='start')
    ax.plot(w, data, '.', label='data')
    ax.invert_xaxis()
    ax.set_title('MPL Fit')
    ax.legend()
    plt.show()
    plt.tight_layout()


#%% vectorization backward
if Z.shape[1] > 1:
    
    V_amide_stack = np.zeros((V_m_z.shape[1]), dtype=float) # we initialize a 2D array filled with zeros and with the same shape as V_m_z
    V_amide_stack[mask_idx] = amide_stack # we assign the values from dB0_stack to B0_stack
    V_amide_stack = V_amide_stack.reshape(V.shape[0], V.shape[1], V.shape[2]) # we reshape B0_stack into a 3D volume
    
    V_rNOE_stack = np.zeros((V_m_z.shape[1]), dtype=float) # we initialize a 2D array filled with zeros and with the same shape as V_m_z
    V_rNOE_stack[mask_idx] = rNOE_stack # we assign the values from dB0_stack to B0_stack
    V_rNOE_stack = V_rNOE_stack.reshape(V.shape[0], V.shape[1], V.shape[2]) # we reshape B0_stack into a 3D volume
   
    V_MT_stack = np.zeros((V_m_z.shape[1]), dtype=float) # we initialize a 2D array filled with zeros and with the same shape as V_m_z
    V_MT_stack[mask_idx] = MT_stack # we assign the values from dB0_stack to B0_stack
    V_MT_stack = V_MT_stack.reshape(V.shape[0], V.shape[1], V.shape[2]) # we reshape B0_stack into a 3D volume
    
    Z_fit_stack = np.zeros((V_m_z.shape[1],np.size(w)), dtype=float)
    Z_fit_stack[mask_idx,:] = FIT_stack
    Z_fit_stack = Z_fit_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(w))
    
    Z_stack = np.zeros((V_m_z.shape[1],np.size(w)), dtype=float)
    Z_stack[mask_idx,:] = Z.transpose()
    Z_stack = Z_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(w))


# %% ==================
# 5) Plot Parametric Maps from Z(3.5 ppm) and MTRasym(3.5ppm)
# =====================
if data_flag == 'real_data':
    slice_of_interest = 2  # pick slice for Evaluation

    plt.figure(figsize=(10, 4))
    plt.subplot(1, 3, 1)
    plt.imshow(V_amide_stack[:, :, slice_of_interest], vmin=0.0, vmax=0.1)
    plt.colorbar()
    plt.title("amide")
    plt.subplot(1, 3, 2)
    plt.imshow(V_rNOE_stack[:, :, slice_of_interest],vmin=0.0,vmax=0.2)
    plt.colorbar()
    plt.title('rNOE')
    plt.subplot(1, 3, 3)
    plt.imshow(V_MT_stack[:, :, slice_of_interest],vmin=0.0,vmax=0.2)
    plt.colorbar()
    plt.title('MT')
    plt.show()
