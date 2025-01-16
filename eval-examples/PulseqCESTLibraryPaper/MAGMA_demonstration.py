import os

# Paths to the target files
current_dir = os.path.dirname(os.path.abspath(__file__))

APTw_file_path = os.path.join(
    current_dir,
    "../seq-library/APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor/APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.py"
)

WASABI_file_path = os.path.join(
    current_dir,
    "../seq-library/WASABI_3T_001_3p7uT_1block_5ms/WASABI_3T_001_3p7uT_1block_5ms.py"
)

# Define functions to load and execute the files on demand
def APTw_seq_gen():
    """Generate the APTw sequence."""
    with open(APTw_file_path, 'r') as f:
        code = f.read()
        exec(code, globals())  # Execute the file in the current global namespace

def WASABI_seq_gen():
    """Generate the WASABI sequence."""
    with open(WASABI_file_path, 'r') as f:
        code = f.read()
        exec(code, globals())  # Execute the file in the current global namespace

# ---- MAGMA 2025 Demonstration ----

# This code file demonstrates the usage of the Pulseq-CEST Library for creating sequence plots, Z-spectra, MTRasym spectra, 
# and parametric maps for different use cases (e.g. varying sequences such as WASABI, APTw and environments like white 
# matter, creatine, and L-arginine, both real and simulated).

# Under each figure header, uncomment the desired code to run the desired subfigure generation.

# Figure 1
# [not applicable -- figure shows screenshot of Pulseq-CEST Library contents]

# Figure 2a, 3a, 4a
# Uncomment to generate the APTw sequence
# APTw_seq_gen()

# Figure 5a
# Uncomment to generate the WASABI sequence
WASABI_seq_gen()

# Figure 3

# Figure 4

# Figure 5

# Figure 6

# Figure S1

# Figure S2
# [not applicable -- figure shows exchange rate estimation over variations in pH level, sourced from cited work]
