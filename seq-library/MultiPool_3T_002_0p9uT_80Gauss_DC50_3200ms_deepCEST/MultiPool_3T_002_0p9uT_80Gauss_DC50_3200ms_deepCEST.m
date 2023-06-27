% MultiPool_3T_002_0p9uT_80Gauss_DC50_3200ms_deepCEST
% Creates a sequence file for a multi pool protocol for 3T as used for the 3T to 7T deepCEST pipeline by Leonie Hunger
% this is the 0.9 µT protocol from  MultiPool_3T_001_0p9uT_80Gauss_DC50_3200ms_GLINT 
% but with a slightly different offsetlist compared to the original publications.
%
% Moritz Zaiss, 2023

% author name for sequence file
author = 'Moritz Zaiss';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% scanner limits
% see pulseq doc for more ino
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;                  % for H [Hz/uT]

%% sequence definitions
% everything in defs gets written as definition in .seq-file
defs.n_pulses      = 80         ; % number of pulses
defs.tp            = 20.48e-3   ; % pulse duration [s]
defs.td            = 20e-3      ; % interpulse delay [s]
defs.Trec          = 0          ; % recovery time [s]
defs.Trec_M0       = 12         ; % recovery time before M0 [s]
defs.M0_offset     = -300       ; % m0 offset [ppm]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.offsets_ppm= [defs.M0_offset  -100. -50. -40. -30. -20. -10. -9.5 -9.  -8.5 -8.  -7.5 -7.  -6.5 -6.  -5.5 -5.  -4.5 -4.  -3.5 -3.  -2.5 -2.  -1.5 -1.  -0.6 -0.4 -0.2  0.   0.2  0.4   0.6  1.   1.5  2.   2.5  3.   3.5  4.   4.5  5.   5.5  6.   6.5  7.   7.5  8.   8.5  9.   9.5 10. 20.  30.  40.  50. 100. ]
defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - ...
    defs.td ;  % saturation time [s]
defs.B0            = 3               ; % B0 [T]
defs.seq_id_string = seqid           ; % unique seq id

%% two different b1 proocols
defs.B1pa        = 0.9;  % mean sat pulse b1 [uT] (the cw power equivalent will be calculated and written to defs below)

%%
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(defs.seq_id_string,'.seq'); % filename

%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = defs.B1pa*gamma_rad*defs.tp; % flip angle of sat pulse
% create pulseq saturation pulse object

satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', defs.tp,'system',seq.sys,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss

[B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,defs.tp,defs.td,1,gamma_hz);
defs.B1rms = B1rms;

%% loop through zspec offsets
offsets_Hz = defs.offsets_ppm*gamma_hz*defs.B0;

% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == defs.M0_offset*gamma_hz*defs.B0
        seq.addBlock(mr.makeDelay(defs.Trec_M0));
    end
    
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase=0;
    for np = 1:defs.n_pulses
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        seq.addBlock(satPulse) % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        if np < defs.n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(defs.td)); % add delay
        end
    end
    
    if defs.spoiling  % spoiling before readout
        seq.addSpoilerGradients()
    end
    
    seq.addPseudoADCBlock(); % readout trigger event
end

%% write definitions
def_fields = fieldnames(defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
saveSaturationPhasePlot(seq_filename);

%% call standard sim
M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');
figure('Name','Z-asym');
plotSimulationResults(M_z,defs.offsets_ppm,defs.M0_offset);
writematrix(M_z', ['M_z_' seq_filename '.txt']);