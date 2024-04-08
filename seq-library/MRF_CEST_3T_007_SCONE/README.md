# Global deep learning optimized pulse sequence for CEST MRF
A spin-lock preparation was added to the original schedule. A pseudo-ADC block instructs the scanner to perform its readout. 

Note that in the original publication the readout flip angle is varied. This can be implemented using the FA provided in the `MRF_CEST_3T_SCONE_007.txt` file (4th column from the left).

## Description
The folder contains a CEST-MRF protocol for brain imaging: 
- B<sub>1</sub> = 0-4 µT
- T<sub>sat</sub> = 0-4 s
- T<sub>R</sub> = 0-4 s 
- Sat. Pulse: Gaussian pulse train (16 ms each).

## Publication
Cohen O, Otazo R. Global deep learning optimization of chemical exchange saturation transfer magnetic resonance fingerprinting acquisition schedule. NMR in Biomedicine. 2023; 36(10):e4954. doi:10.1002/nbm.4954
(https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/10.1002/nbm.4954#)
