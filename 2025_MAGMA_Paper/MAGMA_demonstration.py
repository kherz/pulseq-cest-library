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

# Figure 5a
# subprocess.run(["python3", "../seq-library/WASABI_3T_001_3p7uT_1block_5ms/WASABI_3T_001_3p7uT_1block_5ms.py", "11.5", "15.5"])

# Figure 3

# Figure 4

# Figure 6a, 6b
# subprocess.run(["python3", "WM_comparison.py"])

# Figure S1

# Figure S2
# [not applicable -- figure shows exchange rate estimation over variations in pH level, sourced from cited work]
