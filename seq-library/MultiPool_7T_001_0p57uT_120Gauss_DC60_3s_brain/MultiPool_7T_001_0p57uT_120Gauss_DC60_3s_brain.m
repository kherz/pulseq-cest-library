%% MultiPool_7T_001_0p57uT_120Gauss_DC60_3s_brain
% Creates a sequence file for an multi pool protocol for 7T according to
% https://cest-sources.org/doku.php?id=standard_cest_protocols
% MP_1 (//GUFI TH2//)

%
% Kai Herz 2020
% kai.herz@tuebingen.mpg.de

% author name for sequence file
author = 'Kai Herz';

%% get id of generation file
if strcmp(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = 120              ; % number of pulses
seq_defs.tp            = 15e-3           ; % pulse duration [s]
seq_defs.td            = 10e-3            ; % interpulse delay [s]
seq_defs.Trec          = 0             ; % recovery time [s]
seq_defs.Trec_M0       = 12             ; % recovery time before M0 [s]
seq_defs.M0_offset     = -300           ; % m0 offset [ppm]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [seq_defs.M0_offset -100 -50 -20 -12 -9 -7.25 -6.25 -5.5 -4.7 -4 -3.3 -2.7...
    -2 -1.7 -1.5 -1.1 -0.9 -0.6 -0.4 0 0.4 0.6 0.95 1.1 1.25 1.4 1.55 1.7 1.85 ...
    2 2.15 2.3 2.45 2.6 2.75 2.9 3.05 3.2 3.35 3.5 3.65 3.8 3.95 4.1 4.25 4.4 ...
    4.7 5.25 6.25 8 12 20 50 100 seq_defs.M0_offset ];  % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
    seq_defs.td ;  % saturation time [s]
seq_defs.B0            = 7               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id


%% get info from struct
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
Trec        = seq_defs.Trec;        % recovery time between scans [s]
Trec_M0     = seq_defs.Trec_M0;     % recovery time before m0 scan [s]
tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B0          = seq_defs.B0;          % B0 [T]

%% two different b1 proocols
B1pa        = 0.6;  % mean sat pulse b1 [uT]
%B1pa        = 0.9;  % mean sat pulse b1 [uT]

%%
spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
lims = Get_scanner_limits();

%% create scanner events
% satpulse
gyroRatio_hz  = 42.5764;                  % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;        % [rad/uT]
fa_sat        = B1pa*gyroRatio_rad*tp; % flip angle of sat pulse
% create pulseq saturation pulse object

satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', tp,'system',lims,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss


[B1cwpe,B1cwae,B1cwae_pure,alpha]= calc_power_equivalents(satPulse,tp,td,1,gyroRatio_hz);
seq_defs.B1cwpe = B1cwpe;

% spoilers
spoilRiseTime = 1e-3;
spoilDuration = 4500e-6+ spoilRiseTime; % [s]
% create pulseq gradient object
[gxSpoil, gySpoil, gzSpoil] = Create_spoiler_gradients(lims, spoilDuration, spoilRiseTime);

% pseudo adc, not played out
pseudoADC = mr.makeAdc(1,'Duration', 1e-3);

%% loop through zspec offsets
offsets_Hz = offsets_ppm*gyroRatio_hz*B0;

% init sequence
seq = mr.Sequence();
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == seq_defs.M0_offset*gyroRatio_hz*B0
        
        seq.addBlock(mr.makeDelay(Trec_M0));
    end
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    accumPhase=0;
    for np = 1:n_pulses
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated pahse from previous rf pulse
        seq.addBlock(satPulse) % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        if np < n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(td)); % add delay
        end
    end
    if spoiling % spoiling before readout
        seq.addBlock(gxSpoil,gySpoil,gzSpoil);
    end
    seq.addBlock(pseudoADC); % readout trigger event
end


%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
save_seq_plot(seq_filename);





