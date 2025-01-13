%% make example combined seq file
% Add pulseq version to path in matlab
%% combine files
% Define path for CEST seq file
prep_fn = 'seq\MultiPool_3T_002_0p9uT_80Gauss_DC50_3200ms_deepCEST.seq';

% Define path for readout seq file
ro_fn = 'seq\gre.seq';

% Define path and name for combined seq file
combined_fn = 'seq\MultiPool_3T_002_0p9uT_80Gauss_DC50_3200ms_deepCEST_gre.seq';

combinePrepAndReadoutPulseqFiles(prep_fn,ro_fn,combined_fn);