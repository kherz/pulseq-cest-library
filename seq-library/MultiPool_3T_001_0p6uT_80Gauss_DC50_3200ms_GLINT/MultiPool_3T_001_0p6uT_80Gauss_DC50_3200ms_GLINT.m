% MultiPool_3T_001_0p6uT_80Gauss_DC50_3200ms_GLINT
% Creates a sequence file for a multi pool protocol for 3T according to
% https://cest-sources.org/doku.php?id=standard_cest_protocols
% 
% Note: This protocol is set to the 0.6uT field amplitude (first scan) by default.
% If you would like to generate the 0.9uT field amplitude (second scan), change the B1pa parameter below.
%
% Lukas Kamm, 2021
% lukas.kamm@gmail.com

% author name for sequence file
author = 'Lukas Kamm';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = 80         ; % number of pulses
seq_defs.tp            = 20.48e-3   ; % pulse duration [s]
seq_defs.td            = 20e-3      ; % interpulse delay [s]
seq_defs.Trec          = 0          ; % recovery time [s]
seq_defs.Trec_M0       = 12         ; % recovery time before M0 [s]
seq_defs.M0_offset     = -300       ; % m0 offset [ppm]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [seq_defs.M0_offset -100  -50  -40  -30  -20  -10  -9.5  -9 -8.5 -8 -7.5 -7 -6.5 -6 -5.5 -5 -4.5 -4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 -0.25 0 0.25 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5  8 8.5 9 9.5 10 20 30 40 50 100 ];  % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
    seq_defs.td ;  % saturation time [s]
seq_defs.B0            = 3               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id

%% two different b1 proocols
B1pa        = 0.6;  % mean sat pulse b1 [uT] (the cw power equivalent will be calculated and written to seq_defs below)
%B1pa        = 0.9;  % mean sat pulse b1 [uT]

%%
spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
seq = SequenceSBB(getScannerLimits());

%% create scanner events
% satpulse
gyroRatio_hz  = 42.5764;                  % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;        % [rad/uT]
fa_sat        = B1pa*gyroRatio_rad*seq_defs.tp; % flip angle of sat pulse
% create pulseq saturation pulse object

satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', seq_defs.tp,'system',seq.sys,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss

[B1cwpe,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,seq_defs.tp,seq_defs.td,1,gyroRatio_hz);
seq_defs.B1cwpe = B1cwpe;

%% loop through zspec offsets
offsets_Hz = seq_defs.offsets_ppm*gyroRatio_hz*seq_defs.B0;

% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == seq_defs.M0_offset*gyroRatio_hz*seq_defs.B0
        seq.addBlock(mr.makeDelay(seq_defs.Trec_M0));
    end
    
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase=0;
    for np = 1:seq_defs.n_pulses
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        seq.addBlock(satPulse) % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        if np < seq_defs.n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(seq_defs.td)); % add delay
        end
    end
    
    if spoiling % spoiling before readout
        seq.addSpoilerGradients()
    end
    
    seq.addPseudoADCBlock(); % readout trigger event
end

%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
saveSaturationPhasePlot(seq_filename);

%% call standard sim
M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_001_bmsim.yaml');