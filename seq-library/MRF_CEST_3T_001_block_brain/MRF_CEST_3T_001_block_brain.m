%% MRF_CEST_3T_001_block_brain
% Creates a sequence file for an MRF-CEST protocol according to Figure 1 of:
% doi:10.1016/j.neuroimage.2019.01.034.
%
% Kai Herz 2020
% kai.herz@tuebingen.mpg.de

% author name for sequence file
author = 'Kai Herz';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = [4 2 4 3]       ; % number of pulses
seq_defs.tp            = 200e-3          ; % pulse duration [s]
seq_defs.td            = 10e-3           ; % interpulse delay [s]
seq_defs.Trec          = 2.5             ; % approx [s]
seq_defs.DCsat         = seq_defs.tp/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [4 3 3.5 10]; % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses.*(seq_defs.tp+ seq_defs.td)-seq_defs.td;
seq_defs.B0            = 3               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id
seq_defs.B1pa           = [1.2 0.8 2 3]   ;


%% get info from struct
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
Trec        = seq_defs.Trec;
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B0          = seq_defs.B0;          % B0 [T]
B1          = seq_defs.B1pa;          % mean sat pulse b1 [uT]
spoiling    = 1;                   % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
seq = SequenceSBB(getScannerLimits());

%% create scanner events
% satpulse
gamma_hz  =seq.sys.gamma*10e-6;                  % for H [Hz/uT]
gamma_rad = gamma_hz*2*pi;        % [rad/uT]


%% loop through zspec offsets
offsets_Hz = offsets_ppm*gamma_hz*B0;

meas_id = 1;
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    
    seq.addBlock(mr.makeDelay(Trec)); % recovery time
    fa_sat        = B1(meas_id)*gamma_rad*tp; % flip angle of sat pulse
    satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', tp, 'system', seq.sys);
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase=0;
    for np = 1:n_pulses(meas_id)
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        seq.addBlock(satPulse) % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        
        if np < n_pulses(meas_id) % delay between pulses
            
            seq.addBlock(mr.makeDelay(td)); % add delay
            
        end
    end
    if spoiling % spoiling before readout
        seq.addSpoilerGradients()
    end
    seq.addPseudoADCBlock(); % readout trigger event
    meas_id = meas_id + 1;
end


%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% write sequence
seq.write(seq_filename);

%% plot
saveSaturationPhasePlot(seq_filename);

%% call standard sim
M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');

