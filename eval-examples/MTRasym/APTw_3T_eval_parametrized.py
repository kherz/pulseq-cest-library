from pathlib import Path
import os
import matplotlib.pyplot as plt
import numpy as np
import pydicom
import pypulseq as pp  # Ensures pypulseq is imported
from bmctool.simulate import simulate
from csaps import csaps
from scipy.interpolate import interp1d

def evaluate_APTw_3T(
    data_flag="real_data",
    data_path="",
    bmsim_filename="WM_3T_001_bmsim.yaml",
    seq_filename="APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
    interpolation="spline",
    dB0_shift=-0.02  # Adjusted for phantom
):
    # Define sequence file path (Ensures Paths are Correct)
    seq_name = Path(seq_filename)
    seq_path = Path.cwd().parent / "seq-library" / seq_name.stem / seq_name
    assert seq_path.is_file(), f"Sequence file not found: {seq_path}"

    # Initialize Sequence Object Before Use
    seq = pp.Sequence()
    seq.read(seq_path)
    
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

    # Normalize Data
    M0_idx = np.where(abs(offsets) >= abs(m0_offset))[0]
    if len(M0_idx) > 0:
        M0 = np.mean(m_z[M0_idx, :], 0)
        offsets = np.delete(offsets, M0_idx)
        m_z = np.delete(m_z, M0_idx, axis=0)
        Z = m_z / M0
    else:
        print("m0_offset not found in offset")

    # Helper Function for Polynomial Evaluation
    def ppval(p, x):
        if callable(p):
            return p(x)
        else:
            n = len(p) - 1
            result = np.zeros_like(x)
            for i in range(n, -1, -1):
                result = result * x + p[i]
            return result

    # Smoothing (Interpolation)
    if interpolation == "linear":
        Z_corr = np.zeros_like(Z)
        w = offsets
        dB0_stack = np.zeros(Z.shape[1])

        for ii in range(Z.shape[1]):
            if np.all(np.isfinite(Z[:, ii])):
                f = interp1d(w, Z[:, ii], kind='linear', fill_value='extrapolate')
                w_fine = np.arange(-1, 1.005, 0.005)
                z_fine = f(w_fine)

                min_idx = np.argmin(z_fine)
                dB0_stack[ii] = w_fine[min_idx]

                Z_corr[:, ii] = f(w + dB0_stack[ii] + dB0_shift)

    elif interpolation == "spline":
        Z_corr = np.zeros_like(Z)
        w = offsets
        dB0_stack = np.zeros(Z.shape[1])

        for ii in range(Z.shape[1]):
            if np.all(np.isfinite(Z[:, ii])):
                pp_spline = csaps(w, Z[:, ii], smooth=0.99)
                w_fine = np.arange(-1, 1.005, 0.005)
                z_fine = ppval(pp_spline, w_fine)

                min_idx = np.argmin(z_fine)
                dB0_stack[ii] = w_fine[min_idx]

                Z_corr[:, ii] = ppval(pp_spline, w + dB0_stack[ii] + dB0_shift)

    # Calculate MTRasym-Spectrum
    Z_ref = Z_corr[::-1, :]
    MTRasym = Z_ref - Z_corr

    # Plot Results
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

# Example Usage
if __name__ == "__main__":
    evaluate_APTw_3T()
