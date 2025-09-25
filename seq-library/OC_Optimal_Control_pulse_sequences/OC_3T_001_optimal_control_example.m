%% OC_3T_001_optimal_control_example
% Creates a sequence file for an optimal control pulse protocol at 3T
% Features optimized B1 waveforms for enhanced CEST contrast
%
% Clemens Stilianu 2024
% stilianu@tugraz.at

% author name for sequence file
author = 'Clemens Stilianu';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% scanner limits
% see pulseq doc for more info
seq = SequenceSBB(getScannerLimits());
gamma_hz = seq.sys.gamma*1e-6; % for H [Hz/uT]

%% sequence definitions
% everything in defs gets written as definition in .seq-file
defs.n_pulses      = 10                ; % number of pulses
defs.tp            = 100e-3           ; % pulse duration [s]
defs.td            = 10e-3            ; % interpulse delay [s]
defs.Trec          = 3.5              ; % recovery time [s]
defs.Trec_M0       = 3.5              ; % recovery time before M0 [s]
defs.M0_offset     = -300             ; % m0 offset [ppm]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.offsets_ppm   = [defs.M0_offset -4:0.25:4]; % offset vector [ppm]
defs.num_meas      = numel(defs.offsets_ppm); % number of repetitions
defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - defs.td; % saturation time [s]
defs.FREQ          = 123.2627         ; % Approximately 3 T [MHz]
defs.B0            = defs.FREQ/(gamma_hz); % Calculate B0 [T]
defs.B1pa          = 1.0              ; % target B1 amplitude [uT]
defs.spoiling      = 1                ; % 0=no spoiling, 1=before readout

defs.seq_id_string = seqid; % unique seq id

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename

%% create scanner events
% optimal control saturation pulse
gamma_rad = gamma_hz*2*pi; % [rad/uT]
fa_sat = defs.B1pa*gamma_rad*defs.tp; % flip angle of sat pulse

% Create OC Pulseq saturation pulse object with low-pass filter
satPulse = makeOCPulse(fa_sat, 'system', seq.sys, 'duration', defs.tp, 'useLowPass', true);

% Calculate power equivalent metrics
[B1rms,B1cwae,B1cwae_pure,alpha] = calculatePowerEquivalents(satPulse,defs.tp,defs.td,1,gamma_hz);
defs.B1rms = B1rms;

% Display pulse characteristics
B1_OC = satPulse.signal/gamma_hz; % in micro T
RMS_OC = rms(satPulse.signal/gamma_hz);
fprintf('OC pulse characteristics:\n');
fprintf('  Duration: %.1f ms\n', defs.tp*1000);
fprintf('  B1 RMS: %.2f µT\n', RMS_OC);
fprintf('  B1 peak: %.2f µT\n', max(abs(B1_OC)));
fprintf('  B1 CWAE: %.2f µT\n', B1cwae);

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
    satPulse.freqOffset = currentOffset; % set frequency offset of the pulse
    accumPhase = 0;
    for np = 1:defs.n_pulses
        satPulse.phaseOffset = mod(accumPhase,2*pi); % set accumulated phase from previous rf pulse
        seq.addBlock(satPulse); % add sat pulse
        % calc phase for next rf pulse
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6),2*pi);
        if np < defs.n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(defs.td)); % add delay
        end
    end
    if defs.spoiling % spoiling before readout
        seq.addSpoilerGradients();
    end
    seq.addPseudoADCBlock(); % readout trigger event
end

%% write definitions
def_fields = fieldnames(defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot sequence
fprintf('Sequence file created: %s\n', seq_filename);
saveSaturationPhasePlot(seq_filename);

%% call standard simulation
% Note: Add YAML file when available
M_z = simulate_pulseqcest(seq_filename, 'OC_APT_phantom_sim_test.yaml');

%% plot results (when simulation is available)
plotSimulationResults(M_z, defs.offsets_ppm, defs.M0_offset);
% writematrix(M_z', ['M_z_' seq_filename '.txt']);

fprintf('OC pulse sequence generation completed successfully.\n');

