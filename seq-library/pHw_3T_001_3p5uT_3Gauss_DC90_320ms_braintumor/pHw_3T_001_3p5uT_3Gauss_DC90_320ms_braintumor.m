%% pHw_3T_001_3p5uT_3Gauss_DC90_320ms_braintumor
% Creates a sequence file for an ph weighted protocol according to
% https://doi.org/10.1002/mrm.27204

% This is a multi-slice method with repeating saturation between slices.
% There is an ADC event for each slice, which means, each frequency offset
% contains nSlices (25) ADC events

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
seq_defs.n_pulses      = 3              ; % number of pulses
seq_defs.tp            = 100e-3           ; % pulse duration [s]
seq_defs.td            = 10e-3            ; % interpulse delay [s]
seq_defs.M0_offset     = []           ; % m0 offset [ppm]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [-3.5:0.1:-2.5, -0.3:0.1:0.3, 2.5:0.1:3.5]; % ?3.5 to ?2.5 ppm, ?0.3 to +0.3 ppm, and +2.5 to +3.5 ppm, all in increments of 0.1 ppm
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)+1   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
    seq_defs.td ;  % saturation time [s]
seq_defs.B0            = 3               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id
seq_defs.nSlices       = 25;  % 


%% get info from struct
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B0          = seq_defs.B0;          % B0 [T]
B1peak      = 6;  % mean sat pulse b1 [uT]
spoiling    = 1;  % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
seq = SequenceSBB(getScannerLimits());

%% create scanner events
% satpulse
gamma_hz  =seq.sys.gamma*10e-6;                  % for H [Hz/uT]
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = gamma_rad*tp; % flip angle of sat pulse
% create pulseq saturation pulse object
satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', tp,'system',seq.sys,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
satPulse.signal = (satPulse.signal)./max(satPulse.signal)*B1peak*gamma_hz;
[B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,tp,td,1,gamma_hz);
seq_defs.B1rms = B1rms;

%% loop through zspec offsets
offsets_Hz = offsets_ppm*gamma_hz*B0;

% unsaturated m0
for nSl = 1:seq_defs.nSlices
    seq.addBlock(mr.makeDelay(seq_defs.Tsat));
    seq.addSpoilerGradients()
    seq.addPseudoADCBlock(); % readout trigger event
end

% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    for nSl = 1:seq_defs.nSlices
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
        seq.addSpoilerGradients()
        seq.addPseudoADCBlock(); % readout trigger event
    end
end


%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);


