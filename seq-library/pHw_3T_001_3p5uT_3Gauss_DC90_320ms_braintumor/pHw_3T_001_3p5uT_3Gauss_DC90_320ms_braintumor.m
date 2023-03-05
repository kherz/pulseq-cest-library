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

%% scanner limits
% see pulseq doc for more ino
seq = SequenceSBB(getScannerLimits());
gamma_hz  =seq.sys.gamma*1e-6;                  % for H [Hz/uT]

%% sequence definitions
% everything in defs gets written as definition in .seq-file
defs.n_pulses      = 3              ; % number of pulses
defs.tp            = 100e-3           ; % pulse duration [s]
defs.td            = 10e-3            ; % interpulse delay [s]
defs.M0_offset     = -300           ; % m0 offset [ppm]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.offsets_ppm   = [-3.5:0.1:-2.5, -0.3:0.1:0.3, 2.5:0.1:3.5]; % ?3.5 to ?2.5 ppm, ?0.3 to +0.3 ppm, and +2.5 to +3.5 ppm, all in increments of 0.1 ppm
defs.num_meas      = numel(defs.offsets_ppm)+1   ; % number of repetition
defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - ...
    defs.td ;  % saturation time [s]
defs.FREQ		   = 127.7292;          % Approximately 3 T 
defs.B0            = defs.FREQ/(gamma_hz);  % Calculate B0    
defs.seq_id_string = seqid           ; % unique seq id
defs.nSlices       = 25;  % 

defs.B1peak      = 6;  % mean sat pulse b1 [uT]
defs.spoiling    = 1;  % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename

%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = gamma_rad*defs.tp; % flip angle of sat pulse
% create pulseq saturation pulse object
satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', defs.tp,'system',seq.sys,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
satPulse.signal = (satPulse.signal)./max(satPulse.signal)*defs.B1peak*gamma_hz;
[B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,defs.tp,defs.td,1,gamma_hz);
defs.B1rms = B1rms;

%% loop through zspec offsets
offsets_Hz = defs.offsets_ppm*defs.FREQ;

% unsaturated m0
for nSl = 1:defs.nSlices
    seq.addBlock(mr.makeDelay(defs.Tsat));
    seq.addSpoilerGradients()
    seq.addPseudoADCBlock(); % readout trigger event
end

% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
    for nSl = 1:defs.nSlices
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
        seq.addSpoilerGradients()
        seq.addPseudoADCBlock(); % readout trigger event
    end
end


%% write definitions
def_fields = fieldnames(defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);
%% sim

M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');

writematrix(M_z', ['M_z_' seq_filename '.txt']);
