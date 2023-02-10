%% DGErho_3T_001_4uT_adiaSL_DC100_0.12s_braintumor
% Creates a sequence file for the DGErho protocol from
% https://doi.org/10.1002/mrm.27857
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
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = 1              ; % number of pulses
seq_defs.tp            = 120e-3         ; % pulse duration [s]
seq_defs.Trec          = 4              ; % recovery time [s]
seq_defs.Trec_M0       = 12             ; % recovery time before M0 [s]
seq_defs.M0_offset     = -300           ; % m0 offset [ppm]
seq_defs.DCsat         =               1; % duty cycle
seq_defs.offsets_ppm   = [seq_defs.M0_offset -299 0.6 0.9 1.2 1.5 -299 0.6 0.9 1.2 1.5]; % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm); % number of repetition
seq_defs.Tsat          = seq_defs.tp + 2*12e-3;  % locking + 2 x adiabatic pulses
seq_defs.FREQ		   = 127.7292          % Approximately 3 T  
seq_defs.B0            = seq_defs.FREQ/(gamma_hz);  % Calculate B0   
seq_defs.seq_id_string = seqid           ; % unique seq id

%% get info from struct
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
Trec        = seq_defs.Trec;        % recovery time between scans [s]
Trec_M0     = seq_defs.Trec_M0;     % recovery time before m0 scan [s]
tp          = seq_defs.tp;          % sat pulse duration [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B1pa        = 4;  % mean sat pulse b1 [uT]
spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename


%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = B1pa*gamma_rad*tp; % flip angle of sat pulse
% create pulseq saturation pulse object
satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', tp, 'system', seq.sys);
adia_SL       = makeSLExpPulses(B1pa, seq.sys);
seq_defs.B1rms = B1pa;


%% loop through zspec offsets
offsets_Hz = offsets_ppm*seq_defs.FREQ; % Z spec offsets [Hz]

% loop through offsets and set pulses and delays
pre_sl = [];
post_sl = [];
accumPhase = 0;
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == seq_defs.M0_offset*seq_defs.FREQ
        if Trec_M0 > 0
            seq.addBlock(mr.makeDelay(Trec_M0));
        end
    else
        if Trec > 0
            seq.addBlock(mr.makeDelay(Trec)); % recovery time
        end
    end
    if currentOffset < 0
        pre_sl = adia_SL{find(ismember(adia_SL(:,2), 'pre_neg')),1};
        post_sl = adia_SL{find(ismember(adia_SL(:,2), 'post_neg')),1};
    else
        pre_sl = adia_SL{find(ismember(adia_SL(:,2), 'pre_pos')),1};
        post_sl = adia_SL{find(ismember(adia_SL(:,2), 'post_pos')),1};
    end
    % set frequency
    pre_sl.freqOffset = currentOffset;
    accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(pre_sl.signal)>0))*1e-6),2*pi);
    
    satPulse.phaseOffset = mod(accumPhase,2*pi);
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
    
    post_sl.phaseOffset = mod(accumPhase,2*pi);
    post_sl.freqOffset = currentOffset;
    for np = 1:n_pulses
        seq.addBlock(pre_sl)
        seq.addBlock(satPulse) % add sat pulse
        seq.addBlock(post_sl)
        if np < n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(t_d)); % add delay
        end
    end
    if spoiling % spoiling before readout
        seq.addSpoilerGradients()
    end
    seq.addPseudoADCBlock(); % readout trigger event
    accumPhase = 0;
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
M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');


