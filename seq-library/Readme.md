# seq-library
This folder cointains published protocols for different CEST applications. In each subfolder you can find the .seq-file containing the sequence definition.

In addition, the subfolder contain a short description, plot and matlab and python files used to generate the .seq-file. 

general pulseq-CEST preparation definitions are as follows

- defs.n_pulses      = 36              ; % number of pulses
- defs.tp            = 50e-3           ; % pulse duration [s]
- defs.td            = 5e-3            ; % interpulse delay between pulses[s]
- defs.Trec          = 3.5             ; % recovery time [s]
- defs.Trec_M0       = 3.5             ; % recovery time before M0 [s]
- defs.M0_offset     = -300            ; % m0 offset [ppm]
- defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
- defs.offsets_ppm   = [seq_defs.M0_offset -4:0.25:4]; % offset vector [ppm]
- defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
- defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - ...
                         seq_defs.td ;  % saturation time [s]
- defs.FREQ		   = 127.7292 ;         % Approximately 3 T
- defs.B0            = defs.FREQ/(seq.sys.gamma*1e-6);   % Calculate B0    
- defs.seq_id_string = seqid           ; % unique seq id

- defs.B1pa        = 1.78;  % mean sat pulse B1 [uT]  B1pa*gamma=flip_angle
- defs.B1rms        = 3.7            ; % B1rms, this gets calculated and is typically not used as input.
- defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

- defs.TI            = [10 6 5 4 3 2.5 2 1.5 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];  % inversion time before ADC, for T1 mapping sequences

