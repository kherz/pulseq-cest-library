# ---- MAGMA 2025 Demonstration ----

# This code file demonstrates the usage of the Pulseq-CEST Library for creating sequence plots, Z-spectra, MTRasym spectra, 
# and parametric maps for different use cases (e.g. varying sequences such as WASABI, APTw and environments like white 
# matter, creatine, and L-arginine, both real and simulated). 

# In order to run this file make sure the command line is located one level down from the main Pulseq-CEST Library. 
# Under each figure header, uncomment the desired code to run the desired subfigure generation.

import subprocess

# Figure 1
# [not applicable -- figure shows screenshot of Pulseq-CEST Library contents]

# Figure 2a/3a/4a
# subprocess.run(["python3", "../seq-library/APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor/APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.py", "3", "6"])

# Figure 2b, 2c, 2d, 2f, 2g
subprocess.run(["python3", 
                "../eval-examples/MTRasym/APTw_3T_eval.py", 
                "real_data", 
                "",
                "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                "spline"])
  # Notes: 
  # Enter the path to your DICOM directory: /Users/alexliebeskind/Desktop/pulseq-cest-library/data/PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # --- use equivalent path above ---
  # Are the DICOM Files acquired with the same protocol parameters from the PulseqCEST Library? [y/n]: y
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100

# Figure 2e

# Figure 2h/4b, 2i/4c, 2j/4d
# subprocess.run(["python3", 
#                 "../eval-examples/MTRasym/APTw_3T_eval.py", 
#                 "re_simulation", 
#                 "",
#                 "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
#                 "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
#                 "spline"])

# Figure 3b, 3c, 3d, 3f, 3g

# Figure 3e

# Figure 3h, 3i, 3j

# Figure 4e, 4f, 4g
# subprocess.run(["python3", 
#                 "../eval-examples/MTRasym/APTw_3T_eval.py", 
#                 "re_simulation", 
#                 "",
#                 "phantoms/creatine/creatine_3T_pH6.4_T22C_bmsim.yaml",
#                 "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
#                 "spline"])

# Figure 4h, 4i, 4j
# subprocess.run(["python3", 
#                 "../eval-examples/MTRasym/APTw_3T_eval.py", 
#                 "re_simulation", 
#                 "",
#                 "WM_3T_001_bmsim.yaml",
#                 "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
#                 "spline"])

# Figure 5a
# subprocess.run(["python3", "../seq-library/WASABI_3T_001_3p7uT_1block_5ms/WASABI_3T_001_3p7uT_1block_5ms.py", "11.5", "15.5"])

# Figure 5b, 5c, 5d, 5e
# TO DO

# Figure 5f, 5g
# TO DO

# Figure 6a, 6b
# subprocess.run(["python3", "WM_comparison.py"])

# Figure S1a, S1b
# TO DO

# Figure S1c, S1d
# TO DO

# Figure S1e, S1f
# TO DO

# Figure S1g, S1h
# TO DO

# Figure S1i, S1j
# TO DO

# Figure S1k, S1l
# TO DO

# Figure S1m, S1n
# TO DO

# Figure S2
# [not applicable -- figure shows exchange rate estimation over variations in pH level, sourced from cited work]
