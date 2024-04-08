# Pulse sequence for CEST MRF
A spin-lock preparation was added to the original schedule. A pseudo-ADC block instructs the scanner to perform its readout. 

Note that the original publication used a 100% duty cycle which is not available for all scanner models and vendors. The duty cycle implemented here is 50% (100 ms 'on'). It can be modified as needed using the provided .py files located in this folder.

## Description
The folder contains two LOAS optimized protocols (of length 10 and 40 images) for semisolid MT MRF brain imaging: 
- B<sub>1</sub> = 0-2 µT
- T<sub>sat</sub> = 0-2 s
- T<sub>Rec</sub> = 0-5 s 
- ω = 8-50 ppm

## Publication
Kang B, Kim B, Park H, Heo HY. Learning-based optimization of acquisition schedule for magnetization transfer contrast MR fingerprinting. NMR Biomed. 2022 May;35(5):e4662. doi: 10.1002/nbm.4662. Epub 2021 Dec 22. PMID: 34939236; PMCID: PMC9761585.
