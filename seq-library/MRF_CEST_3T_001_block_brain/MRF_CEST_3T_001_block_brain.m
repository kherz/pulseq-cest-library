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

%% scanner limits
% see pulseq doc for more ino
seq = SequenceSBB(getScannerLimits());
gamma_hz  =seq.sys.gamma*1e-6;                  % for H [Hz/uT]

%% sequence definitions
% everything in defs gets written as definition in .seq-file
defs.n_pulses      = [4 2 4 3]       ; % number of pulses
defs.tp            = 200e-3          ; % pulse duration [s]
defs.td            = 10e-3           ; % interpulse delay [s]
defs.Trec          = 2.5             ; % approx [s]
defs.DCsat         = defs.tp/(defs.tp+defs.td); % duty cycle
defs.offsets_ppm   = [4 3 3.5 10]; % offset vector [ppm]
defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
defs.Tsat          = defs.n_pulses.*(defs.tp+ defs.td)-defs.td;
defs.FREQ		   = 127.7292;         % Approximately 3 T 
defs.B0            = defs.FREQ/(gamma_hz);  % Calculate B0    
defs.seq_id_string = seqid           ; % unique seq id
defs.B1pa           = [1.2 0.8 2 3]   ;% mean sat pulse b1 [uT]
defs.spoiling    = 1;                   % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename


%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]


%% loop through zspec offsets
offsets_Hz = defs.offsets_ppm*defs.FREQ;

meas_id = 1;
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    
    seq.addBlock(mr.makeDelay(defs.Trec)); % recovery time
    fa_sat        = defs.B1pa(meas_id)*gamma_rad*defs.tp; % flip angle of sat pulse
    satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', defs.tp, 'system', seq.sys);
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase=0;
    for np = 1:defs.n_pulses(meas_id)
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        seq.addBlock(satPulse) % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        
        if np < defs.n_pulses(meas_id) % delay between pulses
            
            seq.addBlock(mr.makeDelay(defs.td)); % add delay
            
        end
    end
    if defs.spoiling % spoiling before readout
        seq.addSpoilerGradients()
    end
    seq.addPseudoADCBlock(); % readout trigger event
    meas_id = meas_id + 1;
end


%% write definitions
def_fields = fieldnames(defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% write sequence
seq.write(seq_filename);

%% plot
saveSaturationPhasePlot(seq_filename);

%% call standard sim
M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');
writematrix(M_z', ['M_z_' seq_filename '.txt']);
