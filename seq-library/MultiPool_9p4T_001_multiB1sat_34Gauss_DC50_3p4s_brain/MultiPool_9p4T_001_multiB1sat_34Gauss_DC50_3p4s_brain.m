%% MultiPool_9p4T_001_multiB1sat_34Gauss_DC50_3p4s_brain
% Creates a sequence file for an multi pool protocol for 9.4T according to
% https://cest-sources.org/doku.php?id=standard_cest_protocols
% MP_1 (//GUFI TH2//)

% author name for sequence file
author = 'Felix Glang';

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
seq_defs.n_pulses      = 34              ; % number of pulses
seq_defs.tp            = 50e-3           ; % pulse duration [s]
seq_defs.td            = 50e-3            ; % interpulse delay [s]
seq_defs.Trec          = 3.5             ; % recovery time [s]
seq_defs.Trec_M0       = 12             ; % recovery time before M0 [s]
seq_defs.M0_offset     = -300           ; % m0 offset [ppm]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [seq_defs.M0_offset -50.00 -30.00 -20.00 -12.00 -8.00 -6.00 -5.60...
    -5.30 -5.00 -4.70 -4.40 -4.10 -3.80 -3.50 -3.20 -2.90 -2.60 -2.30 -2.10 -2.00 -1.70 ...
    -1.40 -1.15 -1.00 -0.85 -0.75 -0.50 -0.20 0.00 0.20 0.50 0.60 0.70 0.80 1.00 1.15 1.30...
    1.45 1.60 1.75 1.90 2.00 2.10 2.25 2.40 2.55 2.70 2.85 3.00 3.15 3.30 3.40 3.50 3.60 3.70...
    3.85 4.00 4.15 4.30 4.45 4.60 4.75 4.90 5.00 5.15 5.30 5.45 6.25 8.00 12.00 20.00 30.00 50.00 ];  % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
    seq_defs.td ;  % saturation time [s]
seq_defs.FREQ		   = 400.2182;          % Approximately 9.4 T 
seq_defs.B0            = seq_defs.FREQ/(gamma_hz);  % Calculate B0   
seq_defs.seq_id_string = seqid           ; % unique seq id


%% get info from struct
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
Trec        = seq_defs.Trec;        % recovery time between scans [s]
Trec_M0     = seq_defs.Trec_M0;     % recovery time before m0 scan [s]
tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))

%% multiple saturation pulse amplitudes
B1pa        = [0.6, 0.9, 1.1, 1.4];  % mean sat pulse b1 [uT]

%%
spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]


%% loop through sequence
offsets_Hz = offsets_ppm*seq_defs.FREQ;

for currentB1sat = B1pa % loop through B1sat's

    fa_sat        = currentB1sat*gamma_rad*tp; % flip angle of sat pulse
    % create pulseq saturation pulse object
    satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', tp,'system',seq.sys,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
    [B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,tp,td,1,gamma_hz);
    seq_defs.B1rms = B1rms;
    
    % loop through offsets and set pulses and delays
    for currentOffset = offsets_Hz
        if currentOffset == seq_defs.M0_offset*seq_defs.FREQ
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
        seq.addSpoilerGradients()
    end
    seq.addPseudoADCBlock(); % readout trigger event
    end
    fprintf('a B1sat is done\n')
end

%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
saveSaturationPhasePlot(seq_filename);






