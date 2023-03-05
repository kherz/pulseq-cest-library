% WASABI_7T_001_3p7uT_1block_5ms
% Creates a sequence file for a WASABI protocol with 31 offsets and one M0 image at 3T according to:
% https://doi.org/10.1002/mrm.26133
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
defs.n_pulses      = 1              ; % number of pulses
defs.B1rms        = 3.7            ; % b1 for 1 block is cqpe
defs.tp            = 5e-3           ; % pulse duration [s]
defs.Trec          = 3              ; % recovery time [s]
defs.Trec_M0       = 12             ; % recovery time before M0 [s]
defs.M0_offset     = -300           ; % m0 offset [ppm]
defs.offsets_ppm   = [defs.M0_offset linspace(-1.5, 1.5, 31)]; % offset vector [ppm]
defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
defs.Tsat          = defs.tp     ;  % saturation time [s]
defs.FREQ		   = 298.0348         % Approximately 7 T
defs.B0            = defs.FREQ/(gamma_hz);  % Calculate B0     
defs.seq_id_string = seqid           ; % unique seq id
defs.B1          = defs.B1rms;      % B1 [uT]
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(defs.seq_id_string,'.seq'); % filename

%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = defs.B1*gamma_rad*defs.tp; % flip angle of sat pulse

% create pulseq saturation pulse object
satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', defs.tp, 'system', seq.sys); % block pulse

%% loop through zspec offsets
offsets_Hz = defs.offsets_ppm*defs.FREQ;

% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == defs.M0_offset*defs.FREQ
        if defs.Trec_M0 > 0
            seq.addBlock(mr.makeDelay(defs.Trec_M0));
        end
    else
        if defs.Trec > 0
            seq.addBlock(mr.makeDelay(defs.Trec)); % recovery time
        end
    end
    % add single pulse
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    seq.addBlock(satPulse) % add sat pulse
    if defs.spoiling % spoiling before readout
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

plotSimulationResults(M_z,defs.offsets_ppm, defs.M0_offset);
writematrix(M_z', ['M_z_' seq_filename '.txt']);


