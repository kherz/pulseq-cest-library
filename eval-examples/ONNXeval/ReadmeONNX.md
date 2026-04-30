# ONNX-Based CEST Parameter Estimation

This README describes how to use the offline ONNX parameter estimation pipeline for CEST MRI data.

## Overview

The pipeline performs voxel-wise parameter estimation from CEST Z-spectra using a pre-trained ONNX neural network model. The estimation runs entirely offline without requiring GPU access or a Python deep learning environment at inference time.

## Required Files

To run the evaluation, you need the following files:

- **CEST DICOM data** — the measured CEST images exported in **Enhanced DICOM** format (one multi-frame DICOM file per acquisition, not classic slice-by-slice DICOMs)
- **`.onnx` model file** — the pre-trained neural network for parameter estimation
- **`.ini` configuration file** — contains acquisition-specific settings such as the saturation frequency offsets (ppm list), B0 field strength, and normalization parameters that match the training configuration of the ONNX model

## Repository

The example script is hosted at:

**[CEST / Offline ONNX Parameter Estimation · GitLab (FAU RRZE)](https://gitlab.rrze.fau.de/cest/offline_onnx_parameter_estimation)**

Clone the repository before proceeding:

```bash
git clone https://gitlab.rrze.fau.de/cest/offline_onnx_parameter_estimation.git
cd offline_onnx_parameter_estimation
```

## Usage

1. **Export your CEST measurement** from the scanner in **Enhanced DICOM** format. Make sure the DICOM export includes all dynamics (offsets) as frames within a single multi-frame file.

2. **Obtain the matching `.onnx` and `.ini` files** for your specific CEST protocol. These files must correspond to the same saturation scheme and acquisition parameters used during measurement. Contact your local CEST group if you do not have these files.

3. **Place all required files** in a location accessible to the example script:
   ```
   your_working_directory/
   ├── your_cest_measurement.dcm   # Enhanced DICOM
   ├── model.onnx                  # ONNX model
   └── config.ini                  # Configuration file
   ```

4. **Open the example script** provided in the repository and set the paths to your three input files at the top of the script (DICOM path, ONNX path, INI path).

5. **Run the script.** The output will be a set of parameter maps (e.g., CESTR, linewidth, water shift) saved to the specified output directory.

## Notes

- Only **Enhanced DICOM** is supported. Classic DICOM series (one file per slice) are not compatible with the current reader.
- The `.ini` file and `.onnx` model are protocol-specific. Using a model trained on a different offset list or B0 field will produce incorrect results.
- Refer to the repository's own documentation and example notebooks for detailed parameter descriptions and output format.
