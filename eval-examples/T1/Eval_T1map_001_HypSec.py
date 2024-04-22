# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 16:21:32 2023

@author: kouemoin
"""

import numpy as np
import matplotlib.pyplot as plt
import pypulseq as pp
import os
import  pydicom
from scipy.optimize import curve_fit
from bmctool.simulate import simulate

seq = pp.Sequence()

#%% read in associated seq file from Pulse-CEST library 

seq_fn = 'T1map_001_3HypSec'
seq_path = '../../seq-library/'+seq_fn+'/'+seq_fn+'.seq'
seq.read(seq_path)

#%% extract values from the sequence
TI = seq.get_definition('TI') * 1000  # we define the recovery time(TI) in [ms] using the sequence module get_definition and passing the key as string
offsets =TI
t = TI[1:,]
m0_offset = offsets[0] # we assign to the variable m0_offset the first element of offset
Nmeas = int(seq.get_definition('num_meas')[0]) #number of repetition
#%% 2a) read in data from simulation

txt_path = '../../seq-library/'+seq_fn+'/M_z_'+seq_fn+'.seq.txt' 
m_z = np.loadtxt(txt_path);
m_z = np.expand_dims(m_z, axis=1) # we convert a 1D array into a 2D column vector

# %% 2b) re-simulate using a ymal file
# we assume you are in the path of the present file

config_path = '../../sim-library/WM_3T_default_7pool_bmsim.yaml';
sim = simulate(config_file=config_path, seq_file=seq_path) # we simulate the sequence using the sequence file and yaml file
m_z = sim.get_zspec()[1]
m_z = np.expand_dims(m_z, axis=1) # we convert a 1D array into a 2D column vector

#%% 2c) read data from measurement (dicom)

path = '../example_data/dcm/PULSEQ_HYBRID_GRE_2_2_5_T1/'  # we define a variable path and we assign it a string containing the file path to a directory where DICOM files are located

collection = [pydicom.dcmread(path + filename) for filename in sorted(os.listdir(path))] # we create a collection of DICOM objects from the files in the directory

V = np.stack([dcm.pixel_array for dcm in collection]) # we create an array and we use list comprehension to extract the pixel arrays from each dicom objects in collection
V = np.transpose(V, (1, 2, 0)) # we trnspose V into the right shape
sz = V.shape
V = np.reshape(V, [sz[0], sz[1], Nmeas, sz[2] // Nmeas]).transpose(0,1,3,2) # we reshape V into a 4D volume
#%% vectorization
mask = np.squeeze(V[:,:,:,0]) > 100 # we define a mask by selecting all the elements from the first offset and by removing any single-dimensional entries from the shape of the resulting array, effectively converting it to a 2D array. And we applied a filter.
mask_idx = np.where(mask.ravel())[0]  # we create a 1D aaray mask index
V_m_z =  V.reshape(-1, Nmeas).T # we reshape a 4D volume into a 2D matrix and we transpose the resulting matrix
m_z = V_m_z[:, mask_idx]
#%% Normalization
M0_idx = np.where(offsets == m0_offset)[0] # we initialize M0_idx taking the index of the first element in offsets array.

if len(M0_idx) > 0:
    M0 = np.mean(m_z[M0_idx,:],0) # we compute the mean along the colunm axis
    offsets = np.delete(offsets, M0_idx) # we delete the element at the indice specified by M0_idx from the offsets
    m_z = np.delete(m_z, M0_idx,axis=0) # similar to the previous line
    Z = m_z / M0 #normalization
else:
    print('m0_offset not found in offset')

#%% t1 relaxation function
def t1_fit_func(t, A, B, T1):  #function used for fitting data to an exponential function model 
       
    y = A - B * np.exp(-t / T1) # we calculate the analytic function of the T1 relaxation

    return y

#%% perform the t1 fitting
t1_map = np.zeros(Z.shape[1]) # we initialize an array t1_map of zeros with a length equal to the number of columns in the Z array
t1_fit = np.zeros((Z.shape[1],np.size(t)))
for ii in range(Z.shape[1]):
     if np.all(np.isfinite(Z[:, ii])): # we check if all values in the current column ii of Z are finite (i.e., not NaN or infinity)
         
         p0=([0.5, 1, 200]) # we initialize an initial guess. These values are initial estimates for the parameters A , B and T1 of the exponential function t1_fit_func
         try:
              p = curve_fit(t1_fit_func, t, Z[:,ii],  p0=p0, full_output=True, method='lm') # we use curve_fit function to fit t1_fit_func using Levenberg-Marquardt optimization method.
              print(p[0])
              t1_map[ii]= p[0][2] # we extract the fitted T1 values
              t1_fit[ii,:] = t1_fit_func(t, p[0][0],p[0][1], p[0][2])
         except Exception:
            t1_map[ii]= np.nan # in case an exception is raised during the fitting process, we set that particular column to NaN (Not-a-Number)
            t1_fit[ii,:]= np.nan
#%% vectorization backward
if Z.shape[1] > 1:
    T1_stack = np.zeros(( V_m_z.shape[1]), dtype=float) # we initialize a 2D array filled with zeros and with the same shape as V_m_z.
    T1_stack[mask_idx] = t1_map  # we assign the values from t1_map to T1_stack
    T1_reshaped = T1_stack.reshape(V.shape[0], V.shape[1], V.shape[2]) # we reshape T1_stack into a 4D volume
    
    t1_fit_stack = np.zeros((V_m_z.shape[1],np.size(t)), dtype=float)
    t1_fit_stack[mask_idx,:] = t1_fit
    t1_fit_stack = t1_fit_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(t))
    
    Z_stack = np.zeros((V_m_z.shape[1],np.size(t)), dtype=float)
    Z_stack[mask_idx,:] = Z.transpose()
    Z_stack = Z_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(t))
    
    
#%% display plot

fit=t1_fit_stack[51,63,5,:]
data=Z_stack[51,63,5,:]

plt.figure()
plt.plot(t, data, '.', label='data')
plt.plot(t, fit, label='fit')
plt.xlabel('TI [ms]')
plt.title('T1 fit')
plt.legend()
plt.show()         
            
#%% display the image 
slice_of_interest = 5   #pick slice for Evaluation

plt.figure()
plt.imshow(T1_reshaped[:, :, slice_of_interest], vmin=0, vmax=4000, cmap='gray')
plt.colorbar()
plt.title('T1 map [ms] ')           
plt.show()
