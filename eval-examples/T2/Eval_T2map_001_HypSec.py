# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 15:48:58 2023

@author: kouemoin
"""

import numpy as np
import matplotlib.pyplot as plt
import pypulseq as pp
import os
import  pydicom
from scipy.optimize import curve_fit

seq = pp.Sequence()

#%% read die sequence 
seq_path = 'W:/radiologie/mr-physik-data/Mitarbeiter/kouemo/evalualtion-skripts/T2/T2map_001_T2prep.seq' # we define a variable seq_path and we assign it a string containing the location of the sequence

seq.read(seq_path) 

#%% extract some initial values from the sequence
TE = seq.get_definition('TE') # we define the echo time(TE) in [s] using the sequence module get_definition and passing the key as string
offsets = TE
t = TE[1:,]
m0_offset = offsets[0] # we assign to the variable m0_offset the first element of offset
Nmeas = int(seq.get_definition('num_meas')[0]) # number of repetition



#%% read data from measurement (dicom)
path = 'W:/radiologie/mr-physik-data/Mitarbeiter/kouemo/evalualtion-skripts/T2DicomData' # we define a variable path and we assign it a string containing the file path to a directory where DICOM files are located
os.chdir(path) # we change the current working directory of our Python script to the directory specified in the path variable.

collection = [pydicom.dcmread(filename) for filename in os.listdir(path)] # we create a collection of DICOM objects from the files in the directory
V = np.stack([dcm.pixel_array for dcm in collection]) # we create an array and we use list comprehension to extract the pixel arrays from each dicom objects in collection
V = np.transpose(V, (1, 2, 0)) # we transpose V to have the right shape.
sz = V.shape
V = np.reshape(V, [sz[0], sz[1], Nmeas, sz[2] // Nmeas]).transpose(0,1,3,2) # we reshape V into a 4D volume
#%% vectorization
mask = np.squeeze(V[:,:,:,0]) > 100 # we define a mask by selecting all the elements from the first offset and by removing any single-dimensional entries from the shape of the resulting array, effectively converting it to a 2D array. And we applied a filter.
mask_idx = np.where(mask.ravel())[0] # we create a 1D aaray mask index
V_m_z =  V.reshape(-1, Nmeas).T # we reshape a 4D volume into a 2D matrix and we transpose the resulting matrix.
m_z = V_m_z[:, mask_idx] 
#%% Normalization
M0_idx = np.where(offsets == m0_offset)[0] # we initialize M0_idx taking the index of the first element in offsets array.

if len(M0_idx) > 0:
    M0 = np.mean(m_z[M0_idx,:],0) # we compute the mean along the colunm axis
    offsets = np.delete(offsets, M0_idx) # we delete the element at the indices specified by M0_idx from the offsets
    m_z = np.delete(m_z, M0_idx,axis=0) #similar to the previous line
    Z = m_z / M0 # normalization 
else:
    print('m0_offset not found in offset')

#%% t2 relaxation function
def t2_fit_func(t, A, T2): #function used for fitting data to an exponential decay model 
       
    y = A*np.exp(-t / T2) # we calculate the analytic function of the T2 relaxation

    return y

#%% perform the t2 fitting
t2_map = np.zeros(Z.shape[1]) # we initialize an array t2_map of zeros with a length equal to the number of columns in the Z array
t2_fit = np.zeros((Z.shape[1],np.size(t)))
for ii in range(Z.shape[1]):

     if np.all(np.isfinite(Z[:, ii])): # we check if all values in the current column ii of Z are finite (i.e., not NaN or infinity)
          
         p0=([1, 0.06]) # we initialize an initial guess. These values are initial estimates for the parameters A and T2 of the exponential decay function t2_fit_func
         
         try:
              p = curve_fit(t2_fit_func, t, Z[:,ii], p0=p0, full_output=True, method='lm') # we use curve_fit function to fit t2_fit_func and we use Levenberg-Marquardt optimization method.
              print(p[0])
              t2_map[ii]= p[0][1] * 1000 # we extract the fitted T2 values and we converte there to [ms]  
              t2_fit[ii,:] = t2_fit_func(t, p[0][0],p[0][1])
         except Exception:
            t2_map[ii]= np.nan # in case an exception is raised during the fitting process, we set that particular column to NaN (Not-a-Number)
            t2_fit[ii,:]= np.nan
#%% vectorization backward
if Z.shape[1] > 1:
    T2_stack = np.zeros((V_m_z.shape[1]), dtype=float) # we initialize a 2D array filled with zeros and with the same as V_m_z.
    T2_stack[mask_idx] = t2_map # we assign the values from t2_map to T2_stack 
    T2_reshaped = T2_stack.reshape(V.shape[0], V.shape[1], V.shape[2]) # we reshape T2_stack into a 3D volume
    
    t2_fit_stack = np.zeros((V_m_z.shape[1],np.size(t)), dtype=float)
    t2_fit_stack[mask_idx,:] = t2_fit
    t2_fit_stack = t2_fit_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(t))
    
    Z_stack = np.zeros((V_m_z.shape[1],np.size(t)), dtype=float)
    Z_stack[mask_idx,:] = Z.transpose()
    Z_stack = Z_stack.reshape(V.shape[0], V.shape[1],V.shape[2],np.size(t))
#%% display plot

fit = t2_fit_stack[46,31,5,:]
data = Z_stack[46,31,5,:]

plt.figure()
plt.plot(t, data, '.', label='data')
plt.plot(t, fit, label='fit')
plt.xlabel('[s]')
plt.title('T2 fit')
plt.legend()
plt.show()         
            
#%% display the image 

slice_of_interest = 5   #pick slice for Evaluation

plt.figure()
plt.imshow(T2_reshaped[:, :, slice_of_interest], vmin=0, vmax=100, cmap='jet')  # x and y are inverted
plt.colorbar()
plt.title('T2 map [ms] ')           
plt.show()

