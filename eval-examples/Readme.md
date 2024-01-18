# Introduction to eval-examples

Carrying on the idea of pulseq to provide an open framework for development of MR sequences, we tried to establish minimal and easy-to-use evaluation scripts for specific pulseq-cest sequence files. This directory contains codes for both python and MATLAB.

It is organized as follows:

* `MTRasym` - calculates the MTRasymmetry at 3.5ppm (e.g. "APTw_3T_000_2uT_1block_2s_braintumor.seq")
* `T1` - calculates the T1 time map in ms from the saturation recovery experiment (e.g. "T1map_001_3HypSec.seq")
* `T2` - calculates the T2 time map in ms from different TE (e.g. "T2map_001_T2prep.seq")
* `Wasabi` - calculates the the B1 map and B0 map in ppm by fitting the Wasabi equation (e.g. "WASABI_3T_001_3p7uT_1block_5ms.seq")


## Information to DICOM data

DICOM files have to be stored in one data folder with the sequence name, that means, that all ".dcm" or ".IMA" files from the Wasabi seq-file for all slices and offsets have to available in the same folder.

For newer Siemens software versions with enhanced DICOMs (e.g. XA50, XA60) be sure to export them in the 'interoperability' mode.

