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
% everything in defs gets written as definition in .seq-file
defs.n_pulses      = 34              ; % number of pulses
defs.tp            = 50e-3           ; % pulse duration [s]
defs.td            = 50e-3            ; % interpulse delay [s]
defs.Trec          = 3.5             ; % recovery time [s]
defs.Trec_M0       = 12             ; % recovery time before M0 [s]
defs.M0_offset     = -300           ; % m0 offset [ppm]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.offsets_ppm   = [defs.M0_offset -50.00 -30.00 -20.00 -12.00 -8.00 -6.00 -5.60...
    -5.30 -5.00 -4.70 -4.40 -4.10 -3.80 -3.50 -3.20 -2.90 -2.60 -2.30 -2.10 -2.00 -1.70 ...
    -1.40 -1.15 -1.00 -0.85 -0.75 -0.50 -0.20 0.00 0.20 0.50 0.60 0.70 0.80 1.00 1.15 1.30...
    1.45 1.60 1.75 1.90 2.00 2.10 2.25 2.40 2.55 2.70 2.85 3.00 3.15 3.30 3.40 3.50 3.60 3.70...
    3.85 4.00 4.15 4.30 4.45 4.60 4.75 4.90 5.00 5.15 5.30 5.45 6.25 8.00 12.00 20.00 30.00 50.00 ];  % offset vector [ppm]
defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - ...
    defs.td ;  % saturation time [s]
defs.FREQ		   = 400.2182;          % Approximately 9.4 T 
defs.B0            = defs.FREQ/(gamma_hz);  % Calculate B0   
defs.seq_id_string = seqid           ; % unique seq id

%% multiple saturation pulse amplitudes
defs.B1pa        = [0.6, 0.9, 1.1, 1.4];  % mean sat pulse b1 [uT]

%%
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(defs.seq_id_string,'.seq'); % filename

%% loop through sequence
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
offsets_Hz = defs.offsets_ppm*defs.FREQ;

for currentB1sat = defs.B1pa % loop through B1sat's

    fa_sat        = currentB1sat*gamma_rad*defs.tp; % flip angle of sat pulse
    % create pulseq saturation pulse object
    satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', defs.tp,'system',seq.sys,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
    [B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,defs.tp,defs.td,1,gamma_hz);
    defs.B1rms = B1rms;
    
    % loop through offsets and set pulses and delays
    for currentOffset = offsets_Hz
        if currentOffset == defs.M0_offset*defs.FREQ
            seq.addBlock(mr.makeDelay(defs.Trec_M0));
        end
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
    fprintf('a B1sat is done\n')
end

%% write definitions
def_fields = fieldnames(defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
saveSaturationPhasePlot(seq_filename);

M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');

writematrix(M_z', ['M_z_' seq_filename '.txt']);




