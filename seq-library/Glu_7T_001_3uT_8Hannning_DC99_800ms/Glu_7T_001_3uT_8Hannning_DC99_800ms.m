%% Glu_7T_001_3uT_8Hannning_DC99_800ms_braintumor
% Creates a sequence file for an glutamate weighted protocol according to
% https://doi.org/10.1002/mrm.27362 
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
defs.n_pulses      = 8              ; % number of pulses
defs.tp            = 99.8e-3           ; % pulse duration [s]
defs.td            = 0.2e-3            ; % interpulse delay [s]
defs.Trec          = 15             ; %  every  15  s ???
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.M0_offset     = -100;
defs.offsets_ppm   = [defs.M0_offset -20 -4.2:0.2:-1.8 1.8:0.2:4.2 20 -defs.M0_offset]; % 1.8 to 4.2 ppm with a step size of 0.2 ppm
defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - ...
    defs.td ;  % saturation time [s]
defs.FREQ		   = 298.0348;         % Approximately 7 T  
defs.B0            = defs.FREQ/(gamma_hz);  % Calculate B0   
defs.seq_id_string = seqid           ; % unique seq id

defs.spoiling    = 1;  % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename


%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]

fa_sat = deg2rad(3772); % need to find a hanning pulse with that fa
% create pulseq saturation pulse object
satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', defs.tp,'system',seq.sys); % dummy pulse to get the object
hanning_shape = hanning(numel(satPulse.signal));
satPulse.signal = hanning_shape./trapz(satPulse.t,hanning_shape)*(fa_sat./(2*pi));
[B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,defs.tp,defs.td,1,gamma_hz);
defs.B1rms = B1rms;


%% loop through zspec offsets
offsets_Hz = defs.offsets_ppm*defs.FREQ;

% m0 is unsaturated
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    seq.addBlock(mr.makeDelay(defs.Trec)); % recovery time
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
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
writematrix(M_z', ['M_z_' seq_filename '.txt']);
