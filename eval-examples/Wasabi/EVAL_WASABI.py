# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:45:31 2023

@author: kouemoin
"""
import os
from pathlib import Path

import numpy as np
import pydicom
import matplotlib.pyplot as plt
import pypulseq as pp
from scipy.optimize import curve_fit
from bmctool.simulate import simulate

def EVAL_WASABI(data_flag='simulation',  data_path='', bmsim_filename='WM_3T_default_7pool_bmsim.yaml',
                seq_filename='WASABI_3T_001_3p7uT_1block_5ms.seq'):

    
    # Initializations
    seq_name = Path(seq_filename)
    seq_path = Path.cwd().parent.parent / "seq-library" / seq_name.stem / seq_name
    assert seq_path.is_file(), "seq file not found"
    
    
    # Initializations
    seq = pp.Sequence()
    seq.read(seq_path)
    
    offsets = seq.get_definition('offsets_ppm')  # offset vector [ppm]
    m0_offset = seq.get_definition('M0_offset') # corresponds to the first element in the offset vector
    freq = seq.get_definition('FREQ') # frequency [MHz]
    B1 = seq.get_definition('B1') # excitation field (B1) peak amplitude [µT]
    gamma_ = 42.578 
    t_p = seq.get_definition('tp') # pulse duration [s]
    w = offsets[1:]
    Nmeas = len(offsets) # number of repetition
    
    if data_flag == 'simulation':
        # %% 2a)  read in data from simulation
        txt_path = Path.cwd().parent.parent / "seq-library" / seq_name.stem / f'M_z_{seq_name.stem}.seq.txt'
        assert txt_path.is_file(), "Simulation data file not found"
        m_z = np.loadtxt(txt_path)
        m_z = np.expand_dims(m_z, axis=1) # Convert 1D array to 2D column vector
    
    elif data_flag == 're_simulation':
   
        #  %% 2b) Re-simulate using a yaml file
        config_path = Path.cwd().parent.parent / "sim-library" / bmsim_filename
        assert config_path.is_file(), "Config file not found"
        sim = simulate(config_file=config_path, seq_file=seq_path)
        m_z = sim.get_zspec()[1]
        m_z = np.expand_dims(m_z, axis=1)     
  
    elif data_flag == 'real_data':
        
        # %% 2c) Read the WASABI DICOM file and vectorize the data volume
        if data_path == '':
            dcmpath = input('Enter the path to your DICOM directory: ')
        else:
            dcmpath = data_path

        dcmpath = Path(dcmpath)
        assert dcmpath.is_dir(), "DICOM directory not found"
           
        collection = [pydicom.dcmread(dcmpath / filename) for filename in sorted(os.listdir(dcmpath))]

        V = np.stack([dcm.pixel_array for dcm in collection]) # Extract pixel arrays
        V = np.transpose(V, (2, 1, 0)) # Transpose to correct shape
        sz = V.shape
        V = np.reshape(V, [sz[0], sz[1], Nmeas, sz[2] // Nmeas]).transpose(0, 1, 3, 2) # Reshape to 4D volume
    
        #%% vectorization 
        mask = np.squeeze(V[:,:,:,0]) > 100 # we define a mask by selecting all the elements from the first offset and by removing any single-dimensional entries from the shape of the resulting array, effectively converting it to a 2D array. And we applied a filter.
        mask_idx = np.where(mask.ravel())[0]  # we create a 1D aaray mask index
        V_m_z =  V.reshape(-1, Nmeas).T # we reshape a 4D volume into a 2D matrix and we transpose the resulting matrix
        m_z = V_m_z[:, mask_idx]
    
    #%% normalization
    M0_idx = np.where(offsets == m0_offset)[0] # we initialize M0_idx with the index of the first element in offsets vector
    
    if len(M0_idx) > 0:
        M0 = np.mean(m_z[M0_idx,:],0)  # we compute the mean along the colunm axis
        offsets = np.delete(offsets, M0_idx) # we delete the element at the indice specified by M0_idx from the offsets
        m_z = np.delete(m_z, M0_idx,axis=0)  # similar to the previous line
        Z = m_z / M0 # Normalization
    else:
        print('m0_offset not found in offset')
    
    #%% compute the analytic function of wasabi curve    
    def wasabi_fit_2abs(w,B1,offset,c,d):  #function used for fitting data to an sinusoidal function model 
           
        y = np.abs(c - d * np.sin(np.arctan((B1 / (freq / gamma_)) / (w - offset)))**2 *
                   np.sin(np.sqrt((B1 / (freq / gamma_))**2 + (w - offset)**2) * freq * (2 * np.pi) * t_p / 2)**2) # we compute the analytic function of the wasabi curve
    
        return y
    
    # %% perform the Wasabi fit 
    dB0_stack = np.zeros(Z.shape[1]) # we initialize an array dB0_stack of zeros with a length equal to the number of columns in the Z array
    rB1_stack = np.zeros(Z.shape[1]) # same as the previous line
    Z_fit = np.zeros((Z.shape[1],np.size(w))) # same as the previous line
    
    
    for ii in range(Z.shape[1]):
         if np.all(np.isfinite(Z[:, ii])): # we check if all values in the current column ii of Z are finite (i.e., not NaN or infinity)
             
             p0=([3.7, -0.1, 1, 2]) # we initialize an initial guess. These values are initial estimates for the parameters B1, offset, c and d of wasabi_fit_2abs
             try:
                  p = curve_fit(wasabi_fit_2abs, w, Z[:,ii],  p0=p0, full_output=True, method='lm') # we use curve_fit function to fit wasabi_fit_function using Levenberg-Marquardt optimization method.
                  print(p[0])
                  rB1_stack[ii]= p[0][0] / B1 # we extract rB1_stack and normalize it by dividing with B1
                  dB0_stack[ii]= p[0][1] # we extract dB0_stack 
                  Z_fit[ii,:]=wasabi_fit_2abs(w, p[0][0], p[0][1], p[0][2], p[0][3])
                  
             except Exception:
                rB1_stack[ii]= np.nan  # in case an exception is raised during the fitting process, we set that particular column to NaN (Not-a-Number)
                dB0_stack[ii]= np.nan
                Z_fit[ii,:]= np.nan
    
    # %% vectorization backward
     
    if Z.shape[1] > 1:
        
        B0_stack = np.zeros((V_m_z.shape[1]), dtype=float) # we initialize a 2D array filled with zeros and with the same shape as V_m_z
        B0_stack[mask_idx] = dB0_stack # we assign the values from dB0_stack to B0_stack
        B0_reshaped = B0_stack.reshape(V.shape[0], V.shape[1], V.shape[2]) # we reshape B0_stack into a 3D volume
        
        B1_stack = np.zeros((V_m_z.shape[1]), dtype=float)
        B1_stack[mask_idx] = rB1_stack
        B1_reshaped = B1_stack.reshape(V.shape[0], V.shape[1], V.shape[2])
        
        Z_fit_stack = np.zeros((V_m_z.shape[1],np.size(w)), dtype=float)
        Z_fit_stack[mask_idx,:] = Z_fit
        Z_fit_stack = Z_fit_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(w))
        
        Z_stack = np.zeros((V_m_z.shape[1],np.size(w)), dtype=float)
        Z_stack[mask_idx,:] = Z.transpose()
        Z_stack = Z_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(w))
        
       
    #%% Visualization 
    
    if data_flag == 'simulation' or data_flag == 're_simulation':
        
        plt.figure(figsize=(10, 4))
        plt.subplot(1, 2, 1)
        plt.plot(w, np.mean(Z, axis=1), "r.-")  # Mittelwert über Achse 1
        plt.gca().invert_xaxis()  # x-Achse umkehren
        plt.title("Mean Z-spectrum")
        plt.xlabel("Offsets (ppm)")
        plt.ylabel("Normalized Signal")
        plt.grid(True)
        plt.show()
        
       
    if data_flag == 'real_data':
        slice_of_interest = 5  # Pick slice for evaluation
        plt.figure(figsize=(10, 4))
        plt.subplot(1, 2, 1)
        plt.imshow(B0_reshaped[:, :, slice_of_interest].squeeze(), vmin=-0.2, vmax=0.2)
        plt.colorbar()
        plt.title("B0 map")

        plt.subplot(1, 2, 2)
        plt.imshow(B1_reshaped[:, :, slice_of_interest].squeeze(), vmin=0.8, vmax=1.2)
        plt.colorbar()
        plt.title("B1 map")
        plt.show()


if __name__ == "__main__":
    EVAL_WASABI(data_flag='simulation',  
                data_path='', 
                bmsim_filename='WM_3T_default_7pool_bmsim.yaml',
                seq_filename='WASABI_3T_001_3p7uT_1block_5ms.seq')
    
    
