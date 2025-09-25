# Optimal Control Pulses designed for CEST

Optimal Control pulses for enhanced CEST contrast using `makeOCPulse.m`. similar to use like makeGaussPulse().

These pulses are designed to maximize CEST contrast while beeing robust to B0 inhomogeneities. The pulse can be used flexible with arbitrary B1 level and can be used with pulse duration >= 50 ms and any duty cycle, Tsat, pulse number and offset.

**Files:** `OC_3T_001_optimal_control_example.m` (standard library example), `OC_compare_minimal.m` (OC vs Gaussian comparison)

**Usage:** Run `OC_3T_001_optimal_control_example` or `OC_compare_minimal` after adding your YAML simulation file. Default is `OC_APT_phantom_sim_test.yaml`

### Publications
for further information see Stilianu et al.:

* Generalization of Optimal Control Saturation Pulse Design for Robust and High CEST Contrast

* Enhanced and robust contrast in CEST MRI: Saturation pulse shape design via optimal control

Clemens Stilianu 2024, stilianu@hotmail.com
