from pathlib import Path
import os
import matplotlib.pyplot as plt
import numpy as np
import pydicom
import pypulseq  
from bmctool.simulate import simulate
from csaps import csaps
from scipy.interpolate import interp1d

def evaluate_APTw_3T(
    data_flag="real_data",
    data_path="",
    bmsim_filename="WM_3T_001_bmsim.yaml",
    seq_filename="APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
    interpolation="spline",
    smoothingFactor=0.95,
    dB0_shift=0  # Adjusted for phantom
):
    # Define sequence file path (Ensures Paths are Correct)
    seq_name = Path(seq_filename)
    seq_path = Path.cwd().parent / "seq-library" / seq_name.stem / seq_name
    assert seq_path.is_file(), f"Sequence file not found: {seq_path}"

    # Initialize Sequence Object Before Use
    seq = pypulseq.Sequence()
    seq.read(seq_path)

    ROI = 'n'
    x_min = 0
    x_max = 0
    y_min = 0
    y_max = 0
    
    # Get Sequence Definitions
    m0_offset = seq.get_definition("M0_offset")
    offsets = seq.get_definition("offsets_ppm")
    n_meas = len(offsets)

    if data_flag == 'simulation':
        # Read in data from simulation
        seq_path_base = Path.cwd().parent / "seq-library" / seq_name.stem
        m_z = np.loadtxt(os.path.join(seq_path_base, f'M_z_{seq_name}.txt'))
        m_z = np.expand_dims(m_z, axis=1)

    elif data_flag == 're_simulation':
        # Re-simulate using provided BMSim configuration
        config_path = Path.cwd().parent / "sim-library" / bmsim_filename
        assert config_path.is_file(), f"Config file not found: {config_path}"

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
        # Read DICOM Data
        dcmpath = data_path if data_path else input('Enter the path to your DICOM directory: ')
        os.chdir(dcmpath)

        collection = [pydicom.dcmread(os.path.join(dcmpath, filename)) for filename in sorted(os.listdir(dcmpath))]
        V = np.stack([dcm.pixel_array for dcm in collection])
        V = np.transpose(V, (1, 2, 0))
        sz = V.shape
        V = np.reshape(V, [sz[0], sz[1], n_meas, sz[2] // n_meas]).transpose(0, 1, 3, 2)

        # Vectorization
        mask = np.squeeze(V[:, :, :, 0]) > 100
        mask_idx = np.where(mask.ravel())[0]
        V_m_z = V.reshape(-1, n_meas).T
        m_z = V_m_z[:, mask_idx]

        ROI = input('Do you want to specify an ROI [y/n]: ')
        if ROI == 'y':
            slice_of_interest = int(input('Slice of interest: ')) # 5
            x_min = int(input('ROI x-minimum: ')) # 145
            x_max = int(input('ROI x-maximum: ')) # 150
            y_min = int(input('ROI y-minimum: ')) # 95
            y_max = int(input('ROI y-maximum: ')) # 100

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

    # Normalize Data
    M0_idx = np.where(abs(offsets) >= abs(m0_offset))[0]
    if len(M0_idx) >= 0:
        M0 = np.mean(m_z[M0_idx, :], 0)
        offsets = np.delete(offsets, M0_idx)
        m_z = np.delete(m_z, M0_idx, axis=0)
        Z = m_z / M0  # Normalization
    else:
        print("m0_offset not found in offset")

    # Helper function to evaluate piecewise polynomial
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

    # Smoothing 
    if interpolation == "linear":
       
        # Perform linear interpolation
        Z_corr = np.zeros_like(Z)
        w = offsets
        dB0_stack = np.zeros(Z.shape[1])

        for ii in range(Z.shape[1]):
            if np.all(np.isfinite(Z[:, ii])):
                # Create linear interpolation function
                f = interp1d(w, Z[:, ii], kind='linear', fill_value='extrapolate')
                
                # Interpolate values at fine grid points
                w_fine = np.arange(-1, 1.001, 0.001)
                z_fine = f(w_fine)

                # Find index of minimum value
                min_idx = np.argmin(z_fine)
                dB0_stack[ii] = w_fine[min_idx]

                #dB0_stack = dB0_stack - 0.028

                # Interpolate corrected values
                Z_corr[:, ii] = f(w + dB0_stack[ii] + dB0_shift)

    elif interpolation == "spline": 
        # Perform the smoothing spline interpolation
        Z_corr = np.zeros_like(Z)
        w = offsets
        dB0_stack = np.zeros(Z.shape[1])

        for ii in range(Z.shape[1]):
            if np.all(np.isfinite(Z[:, ii])):
                pp = csaps(w, Z[:, ii], smooth=smoothingFactor)
                w_fine = np.arange(-1, 1.005, 0.005)
                z_fine = ppval(pp, w_fine)

                min_idx = np.argmin(z_fine)
                dB0_stack[ii] = w_fine[min_idx]

                Z_corr[:, ii] = ppval(pp, w + dB0_stack[ii] + dB0_shift)

    # Calc of MTRasym-Spectrum
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

        # ROI-specific smoothing to account for smaller selection and associated spectra
        if interpolation == "linear":
            # Perform linear interpolation
            f = interp1d(w, Z_spectrum, kind='linear', fill_value='extrapolate')
            w_fine = np.arange(-1, 1.001, 0.001)
            z_fine = f(w_fine)
            min_idx = np.argmin(z_fine)
            dB0 = w_fine[min_idx]
            Z_spectrum_corr = f(w + dB0)

            # Correct MTR asymmetry spectrum using the same dB0 shift
            f_mtr = interp1d(w, V_MTRasym_spectrum, kind='linear', fill_value='extrapolate')
            V_MTRasym_spectrum_corr = f_mtr(w + dB0 + dB0_shift)

        elif interpolation == "spline":
            # Perform the smoothing spline interpolation
            pp = csaps(w, Z_spectrum, smooth=0.99)
            w_fine = np.arange(-1, 1.005, 0.005)
            z_fine = ppval(pp, w_fine)
            min_idx = np.argmin(z_fine)
            dB0 = w_fine[min_idx]

            Z_spectrum_corr = ppval(pp, w + dB0)

            # Correct MTR asymmetry spectrum using the same dB0 shift
            pp_mtr = csaps(w, V_MTRasym_spectrum, smooth=0.99)
            V_MTRasym_spectrum_corr = ppval(pp_mtr, w + dB0 + dB0_shift)
        
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

# Example Usage
if __name__ == "__main__":
    evaluate_APTw_3T()
