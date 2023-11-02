 # -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 13:12:35 2023

@author: kouemoin

"""
#requires
#pip install pydicom
#pip install csaps

import numpy as np
import matplotlib.pyplot as plt
import pypulseq as pp
from csaps import csaps
import os
import  pydicom
from bmctool.simulate import simulate

seq = pp.Sequence()
# we assume you are in the path of the present file
# cd('pulseq-cest-library\eval-examples\MTRasym')

# %% 1)  read in associated seq file from Pulse-CEST library
seq_fn ='APTw_3T_000_2uT_1block_2s_braintumor'
seq_path = '../../seq-library/'+seq_fn+'/'+seq_fn+'.seq'  # can be a str or a Path

seq.read(seq_path)

m0_offset = seq.get_definition('M0_offset')

offsets = seq.get_definition('offsets_ppm')

Nmeas = len(offsets)

# %% 2a)  read in data from simulation
txt_path = '../../seq-library/'+seq_fn+'/M_z_'+seq_fn+'.seq.txt'  
m_z = np.loadtxt(txt_path);
m_z = np.expand_dims(m_z, axis=1)

# %% 2b) re-simulate using a ymal file

# we assume you are in the path of the present file
# cd('pulseq-cest-library/sim-library')
config_path ='../../sim-library/WM_3T_default_7pool_bmsim.yaml';
sim = simulate(config_file=config_path, seq_file=seq_path)
m_z = sim.get_zspec()[1]
m_z = np.expand_dims(m_z, axis=1)

# %% 2c)  read data from measurement (dicom)
# move to the directory containing the Dicom files
dcmpath = '../example_data/dcm/PULSEQ_HYBRID_GRE_2_2_5_APTW_001/'

# Question if seq file is still the same
question = input('Are the DICOM Files acquired with the same Protocol Parameters from the PulseqCEST Library? [y/n] ')
if question.lower() == 'n':
    seqfile = input('Please enter the path to your seq file: ')
    seq.read(seqfile)
    offsets = seq.get_definition('offsets_ppm')
    Nmeas = len(offsets)

#read data from dicom directory
collection = [pydicom.dcmread(dcmpath+filename) for filename in os.listdir(dcmpath)]

#extract the volume data
V = np.stack([dcm.pixel_array for dcm in collection])
V = np.transpose(V, (1, 2, 0))
sz = V.shape
V = np.reshape(V, [sz[0], sz[1], Nmeas, sz[2] // Nmeas]).transpose(0,1,3,2)

#Vectorization
mask = np.squeeze(V[:,:,:,0]) > 100
mask_idx = np.where(mask.ravel())[0]
V_m_z =  V.reshape(-1, Nmeas).T
m_z = V_m_z[:, mask_idx]


# %% 3)  Evaluation

M0_idx = np.where(offsets == m0_offset)[0]

if len(M0_idx) > 0:
    M0 = np.mean(m_z[M0_idx,:],0)
    offsets = np.delete(offsets, M0_idx)
    m_z = np.delete(m_z, M0_idx,axis=0)
    Z = m_z / M0 #Normalization
else:
    print('m0_offset not found in offset')


#helper function to evaluate piecewise polinomial 
def ppval(p, x):
    if callable(p):
        return p(x)
    else:
        n = len(p) - 1
        result = np.zeros_like(x)
        for i in range(n, -1, -1):
            result = result * x + p[i]
        return result

#perform the smoothing spline interpolation
Z_corr = np.zeros_like(Z)
w = offsets
dB0_stack = np.zeros(Z.shape[1])
for ii in range(Z.shape[1]):
     if np.all(np.isfinite(Z[:, ii])):
         pp = csaps(w, Z[:,ii], smooth=0.95)
         w_fine = np.arange(-1, 1.005, 0.005)
         z_fine = ppval(pp, w_fine)
            
         min_idx = np.argmin(z_fine)
         dB0_stack[ii] = w_fine[min_idx]
            
         Z_corr[:,ii] = ppval(pp, w + dB0_stack[ii])
    
    
# #calc of MTRasym-Spectrum
Z_ref = Z_corr[::-1,:]
MTRasym = Z_ref - Z_corr

#Vectorization Backwards
if Z.shape[1] > 1:
    
    V_MTRasym = np.zeros((V_m_z.shape[0], V_m_z.shape[1]), dtype=float)
    V_MTRasym[1:, mask_idx] = MTRasym
    V_MTRasym_reshaped = V_MTRasym.reshape(V.shape[3], V.shape[0], V.shape[1], V.shape[2]).transpose(1,2,3,0)
    
    V_Z_corr = np.zeros((V_m_z.shape[0], V_m_z.shape[1]), dtype=float)
    V_Z_corr[1:, mask_idx] = Z_corr
    V_Z_corr_reshaped = V_Z_corr.reshape(V.shape[3], V.shape[0], V.shape[1], V.shape[2]).transpose(1,2,3,0)
    
# %% 4)  Plots and further graphics
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(w, np.mean(Z_corr, axis=1), 'r.-')
plt.gca().invert_xaxis()
plt.title('Mean Z-spectrum')

plt.subplot(1, 2, 2)
plt.plot(w, np.mean(MTRasym, axis=1), 'b.-')
plt.xlim([0, 4])
plt.gca().invert_xaxis()
plt.title('Mean MTRasym-spectrum')
plt.show()

# %% 5) display the images

slice_of_interest = 5   #pick slice for Evaluation
desired_offset = 3.5    
offset_of_interest =np.where(offsets == desired_offset)[0]  #pick offset for Evaluation
w_offset_of_interest = w[offset_of_interest]

plt.figure(figsize=(10, 4))

plt.subplot(1, 2, 1)
plt.imshow(V_Z_corr_reshaped[:, :, slice_of_interest, offset_of_interest], vmin=0.5, vmax=1)
plt.colorbar()
plt.title('Z(Δω) = %.2f ppm' % w_offset_of_interest)

plt.subplot(1, 2, 2)
plt.imshow(V_MTRasym_reshaped[:, :, slice_of_interest, offset_of_interest], vmin=-0.05, vmax=0.05)
plt.colorbar()
plt.title('MTRasym(Δω) = %.2f ppm' % w_offset_of_interest)

plt.show()
 
    
    