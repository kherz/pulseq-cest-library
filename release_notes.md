# pulseq-cest-library release notes

# version before 2023_02-mini-hackathon
mrm-paper version + bugfixes

# version after 2023_02-mini-hackathon
solves issues
#30, #32, #33, #34, #36, #39

- changed all B1cpwe variables to B1rms.
- changed gamma definition- now not defined in the code but taken from pulseq
- also changed "cwpe" to "rms" in all non ".m" files, e.g. python and r…
- Changed normalization scan from -1560 to -300 ppm
- Defining B0 again from the Frequency and gamma, because its needed …
- renamed  seq_defs to defs to make files shorter
- simplified all files 
- added test_all 
- M_z_simulation txt files were generated with pulseq-CEST version: https://github.com/kherz/pulseq-cest/commit/a80f833f5a7521e4197f8487fef073a6286cbc1b



