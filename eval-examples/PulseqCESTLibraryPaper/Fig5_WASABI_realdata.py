# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 13:45:31 2023

@author: kouemoin
"""


import numpy as np
import pydicom
import os
import matplotlib.pyplot as plt
import pypulseq as pp
from scipy.optimize import curve_fit
from bmctool.simulate import simulate
from csaps import csaps

seq = pp.Sequence()

#%% read in associated seq file from Pulse-CEST library

seq_fn = 'WASABI_3T_001_3p7uT_1block_5ms'
seq_path = '/Users/alexliebeskind/Desktop/pulseq-cest-library/seq-library/'+seq_fn+'/'+seq_fn+'.seq'  # can be a str or a Path

seq.read(seq_path)

#%%

offsets = seq.get_definition('offsets_ppm')  # offset vector [ppm]
m0_offset = seq.get_definition('M0_offset') # corresponds to the first element in the offset vector
freq = seq.get_definition('FREQ') # frequency [MHz]
B1 = seq.get_definition('B1') # excitation field (B1) peak amplitude [µT]
gamma_ = 42.578 
t_p = seq.get_definition('tp') # pulse duration [s]
w = offsets[1:]
Nmeas = len(offsets) # number of repetition
ROI = 'n'
x_min = 0
x_max = 0
y_min = 0
y_max = 0

#%% 2a)  read in data from simulation
# txt_path = '/Users/alexliebeskind/Desktop/pulseq-cest-library/seq-library/'+seq_fn+'/M_z_'+seq_fn+'.seq.txt' 
# m_z = np.loadtxt(txt_path);
# m_z = np.expand_dims(m_z, axis=1) # we convert a 1D array into a 2D column vector

# #%% 2b) re-simulate using a ymal file
# #we assume you are in the path of the present file

config_path = '/Users/alexliebeskind/Desktop/pulseq-cest-library/sim-library/phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml';
sim = simulate(config_file=config_path, seq_file=seq_path) # we simulate the sequence using the sequence file and yaml file
m_z = sim.get_zspec()[1]
m_z = np.expand_dims(m_z, axis=1) # we convert a 1D array into a 2D column vector
# Plot m_z
plt.figure(figsize=(5, 4))
plt.plot(w, np.mean(m_z[1:], axis=1), '.-')
plt.ylim(0,1)
plt.xlim(-2,2)
plt.gca().invert_xaxis()
plt.title("Mean Z-spectrum")
plt.show()

#%% 2c) read the wasabi dicom file, create a collection and vectorize the data volume

wsbpath = '/Users/alexliebeskind/Desktop/pulseq-cest-library/data/PULSEQ_HYBRID_GRE_2_2_5_WASABI_RR_0016/'  # we define a variable wsbpath and we assign it a string containing the file path to a directory

collection = [pydicom.dcmread(wsbpath+filename) for filename in sorted(os.listdir(wsbpath))] # we create a collection of DICOM objects from the files in the directory

V = np.stack([dcm.pixel_array for dcm in collection]) # we create an array using list comprehension to extract the pixel arrays from each dicom objects in collection
V = np.transpose(V, (1, 2, 0)) # we transpose V into a right shape
sz = V.shape
V = np.reshape(V, [sz[0], sz[1], Nmeas, sz[2] // Nmeas]).transpose(0,1,3,2) # we reshape V into a 4D volume

#%% vectorization 
mask = np.squeeze(V[:,:,:,0]) > 100 # we define a mask by selecting all the elements from the first offset and by removing any single-dimensional entries from the shape of the resulting array, effectively converting it to a 2D array. And we applied a filter.
mask_idx = np.where(mask.ravel())[0]  # we create a 1D aaray mask index
V_m_z =  V.reshape(-1, Nmeas).T # we reshape a 4D volume into a 2D matrix and we transpose the resulting matrix
m_z = V_m_z[:, mask_idx]

# Plot m_z
plt.figure(figsize=(5, 4))
plt.plot(w, np.mean(m_z[1:], axis=1), ".-")
plt.ylim([0, 1000])
plt.xlim([2, -2])
plt.title("Mean Z-spectrum")
plt.show()

ROI = input('Do you want to specify an ROI [y/n]: ')
if ROI == 'y':
    # x_min = int(input('Upper left x-coordinate: '))
    # x_max = int(input('Bottom right x-coordinate: '))
    # y_min = int(input('Upper left y-coordinate: '))
    # y_max = int(input('Bottom right y-coordinate: '))
    slice_of_interest = 5
    # x_min = 96
    # x_max = 97
    # y_min = 76
    # y_max = 77

    x_min = 145
    x_max = 150
    y_min = 95
    y_max = 100

    # mask_ROI = np.zeros_like(mask)
    # mask_ROI[x_min:x_max, y_min:y_max, slice_of_interest] = True
    # mask_idx = np.where(mask.ravel())[0]
    # V_m_z = V.reshape(-1, Nmeas).T
    # m_z = V_m_z[:, mask_idx]

    # Create mask for the specified ROI
    mask_ROI = np.zeros_like(mask)
    mask_ROI[x_min:x_max, y_min:y_max, slice_of_interest] = True
    mask_idx_ROI = np.where(mask_ROI.ravel())[0]
    V_m_z_ROI = V.reshape(-1, Nmeas).T
    m_z_ROI = V_m_z_ROI[:, mask_idx_ROI]

    # mask_idx = mask_idx_ROI
    # V_m_z = V_m_z_ROI
    # m_z = m_z_ROI
    # # Plot m_z
    # plt.figure(figsize=(5, 4))
    # plt.plot(np.mean(m_z_ROI, axis=1), "r.-")
    # plt.title("Mean Z-spectrum for ROI")
    # plt.show()

    # Plot m_z
    # mean_m_z = np.mean(m_z_ROI, axis=1)
    # std_m_z = np.std(m_z_ROI, axis=1)
    # plt.figure(figsize=(5, 4))
    # plt.errorbar(np.arange(len(mean_m_z)), mean_m_z, yerr=std_m_z, fmt="r.-")
    # plt.title("Mean Z-spectrum ROI")
    # plt.show()
    # plt.figure(figsize=(5, 4))
    # plt.plot(np.mean(m_z_ROI, axis=1), ".-")
    # plt.xlabel(r'$\Delta\omega$ [offset]')
    # plt.ylabel(r'Z($\Delta\omega$)')
    # plt.title("Z-spectrum")
    # plt.show()


#%% normalization
M0_idx = np.where(offsets == m0_offset)[0] # we initialize M0_idx with the index of the first element in offsets vector

if len(M0_idx) > 0:
    M0 = np.mean(m_z[M0_idx,:],0)  # we compute the mean along the colunm axis
    offsets = np.delete(offsets, M0_idx) # we delete the element at the indice specified by M0_idx from the offsets
    m_z = np.delete(m_z, M0_idx,axis=0)  # similar to the previous line
    Z = m_z / M0 # Normalization
else:
    print('m0_offset not found in offset')


#-------------------- (uncomment below to plot L-arginin)
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

        Z_corr[:, ii] = ppval(pp, w + dB0_stack[ii]-.07)


# calc of MTRasym-Spectrum
Z_ref = Z_corr[::-1, :]
MTRasym = Z_ref - Z_corr


# %% ==========================
# 4) Plot MEAN ZSpec and MTRasym from Phantom
# =============================

if ROI == 'n':
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.plot(w, np.mean(Z_corr, axis=1), '.-')
    plt.gca().invert_xaxis()
    plt.ylim(0,1)
    plt.xlim(-2,2)
    plt.gca().invert_xaxis()
    plt.title("Mean Z-spectrum")

    plt.subplot(1, 2, 2)
    plt.plot(w, np.mean(MTRasym, axis=1), '.-')
    #plt.xlim([0, 4])
    plt.gca().invert_xaxis()
    plt.title("Mean MTRasym-spectrum")
    plt.show()
# # --------------------
    
#% compute the analytic function of wasabi curve    
def wasabi_fit_2abs(w,B1,offset,c,d):  #function used for fitting data to an sinusoidal function model 
       
    y = np.abs(c - d * np.sin(np.arctan((B1 / (freq / gamma_)) / (w - offset)))**2 *
               np.sin(np.sqrt((B1 / (freq / gamma_))**2 + (w - offset)**2) * freq * (2 * np.pi) * t_p / 2)**2) # we compute the analytic function of the wasabi curve

    return y

#%%perform the wasabi fit 
dB0_stack = np.zeros(Z.shape[1]) # we initialize an array dB0_stack of zeros with a length equal to the number of columns in the Z array
rB1_stack = np.zeros(Z.shape[1]) # same as the previous line
Z_fit = np.zeros((Z.shape[1],np.size(w))) # same as the previous line

for ii in range(Z.shape[1]):
     if np.all(np.isfinite(Z[:, ii])): # we check if all values in the current column ii of Z are finite (i.e., not NaN or infinity)
         
         p0=([3.7, -0.1, 1, 2]) # we initialize an initial guess. These values are initial estimates for the parameters B1, offset, c and d of wasabi_fit_2abs
         try:
              p = curve_fit(wasabi_fit_2abs, w, Z[:,ii],  p0=p0, full_output=True, method='lm') # we use curve_fit function to fit wasabi_fit_function using Levenberg-Marquardt optimization method.
              # print(p[0])
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
   
#%% display the plots 

fig, ax = plt.subplots()

# with ROI
# fit=Z_fit.transpose()[:,:]
# data=Z[:,:]

# specify ROI
fit=Z_fit_stack[145:150,95:100,5,:]
data=Z_stack[145:150,95:100,5,:]
B0=B0_reshaped[145:150,95:100,5]
B1=B1_reshaped[145:150,95:100,5]
fit = fit.mean(axis=(0, 1))
data = data.mean(axis=(0, 1))
B0 = B0.mean(axis=(0, 1))
B1 = B1.mean(axis=(0, 1))
# B1_plot = np.full(w.shape, B1)
# B0_plot = np.full(w.shape, B0)

print("BO")
print(B0)
print("B1")
print(B1)


#ax.plot(w, B1_plot, label='B1')
ax.plot(w, fit, label='fit')
ax.plot(w, data, '.', label='data')

ax.invert_xaxis()
ax.set_title('WASABI Fit')
ax.set_ylim(0,1)
ax.set_xlim([2,-2])
ax.legend()
plt.show()
plt.tight_layout()
#%% display the image

# without ROI 
slice_of_interest = 5   #pick slice for Evaluation
offset_of_interest = 30 #pick offset for Evaluation
w_offset_of_interest = w[offset_of_interest]

plt.figure(figsize=(10, 4))
ax1 = plt.subplot(1, 2, 1)
plt.imshow(B0_reshaped[:, :, slice_of_interest].squeeze(), vmin=-0.2, vmax=0.2)
plt.colorbar()
plt.title('WASABI B0 map')

x_max, y_max = y_max, x_max
x_min, y_min = y_min, x_min

if ROI == 'y': 
    rect1 = plt.Rectangle((x_min, y_min), (x_max - x_min), (y_max - y_min), linewidth=1, edgecolor='red', facecolor='none')
    ax1.add_patch(rect1)

ax2 = plt.subplot(1, 2, 2)
plt.imshow(B1_reshaped[:, :, slice_of_interest].squeeze(), vmin=0.8, vmax=1.2)
plt.colorbar()
plt.title('WASABI B1 map')

if ROI == 'y': 
    rect2 = plt.Rectangle((x_min, y_min), (x_max - x_min), (y_max - y_min), linewidth=1, edgecolor='red', facecolor='none')
    ax2.add_patch(rect2)

plt.show()

