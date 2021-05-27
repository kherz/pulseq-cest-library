%% T1map_001_3HypSec
% Creates a sequence file for a T1map protocol
%
% Patrick Schuenke 2021
% patrick.schuenke@ptb.de

% author name for sequence file
author = 'Patrick Schuenke';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval') 
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);  
else
    [~, seqid] = fileparts(which(mfilename));
end

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = 3              ; % number of pulses
seq_defs.tp            = 8e-3         ; % pulse duration [s]
seq_defs.Trec          = 1              ; % recovery time [s]
seq_defs.TI            = [10 6 5 4 3 2.5 2 1.5 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];
seq_defs.offsets_ppm   = zeros(1, numel(seq_defs.TI));
seq_defs.num_meas      = numel(seq_defs.offsets_ppm); % number of repetition
seq_defs.seq_id_string = seqid           ; % unique seq id

%% get info from struct
TI          = seq_defs.TI;          % "inversion" times [s]
Trec        = seq_defs.Trec;        % recovery time between scans [s]
tp          = seq_defs.tp;          % sat pulse duration [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B1pa        = 20;  % mean sat pulse b1 [uT]
spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
lims = Get_scanner_limits();

%% create scanner events
% satpulse
gyroRatio_hz  = 42.5764;                  % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;        % [rad/uT]

% create pulseq object
hs_pulse = Generate_HS_HP_pulseq_pulses(20,lims);
seq_defs.B1cwpe = B1pa;

% spoilers
rampTime = 1e-3;
spoilDuration = 4500e-6 + rampTime; % in s

spoilAmplitude = 0.8 .* lims.maxGrad; % in Hz/m % FG: hard coded to 80% of maximum gradient, more should be possible
gxSpoil1=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gySpoil1=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gzSpoil1=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);

spoilAmplitude = -0.7 .* lims.maxGrad;
gxSpoil2=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gySpoil2=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gzSpoil2=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);

spoilAmplitude = 0.6 .* lims.maxGrad;
gxSpoil3=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gySpoil3=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gzSpoil3=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);

% pseudo adc, not played out
pseudoADC = mr.makeAdc(1,'Duration', 1e-3);


%% loop through "inversion" times (TI)
% init sequence
seq = mr.Sequence();

hs_pulse.freqOffset = 0;

for t_prep = TI
    if Trec > 0
        seq.addBlock(mr.makeDelay(Trec)); % recovery time
    end
        
    for np = 1:n_pulses
        seq.addBlock(hs_pulse);
        if mod(np-1,3) == 0
            seq.addBlock(gxSpoil1,gySpoil2,gzSpoil3);
        elseif mod(np-1,2) == 0
            seq.addBlock(gxSpoil3,gySpoil1,gzSpoil2);
        else
            seq.addBlock(gxSpoil2,gySpoil3,gzSpoil1);
        end
    end
    
    % add variable "inversion" time delay
    seq.addBlock(mr.makeDelay(t_prep));
    
    % add pseudo adc / readout trigger event
    seq.addBlock(pseudoADC);
end

%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
save_seq_plot(seq_filename);

%% call standard sim
M_z = Run_pulseq_cest_Simulation(seq_filename,'../../sim-library/GM_3T_001_bmsim.yaml');



