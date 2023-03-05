%% APTw_3T_002_2uT_20SincGauss_DC50_2s_braintumor
% Creates a sequence file for an APTw protocol with Sinc-Gaussian pulses, 50% DC and tsat of 2 s
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
defs.n_pulses      = 20              ; % number of pulses
defs.tp            = 50e-3           ; % pulse duration [s]
defs.td            = 50e-3           ; % interpulse delay [s]
defs.Trec          = 3.5             ; % recovery time [s]
defs.Trec_M0       = 3.5             ; % recovery time before M0 [s]
defs.M0_offset     = -300           ; % m0 offset [ppm]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.offsets_ppm   = [defs.M0_offset -4:0.25:4]; % offset vector [ppm]
defs.num_meas      = numel(defs.offsets_ppm)   ; % number of repetition
defs.Tsat          = defs.n_pulses*(defs.tp+defs.td)- ...
                         defs.td ;  % saturation time [s]
defs.FREQ		   = 127.7292;          % Approximately 3 T
defs.B0            = defs.FREQ/(gamma_hz);   % Calculate B0  
defs.seq_id_string = seqid           ; % unique seq id

defs.B1pa        = 2.4;  % mean sat pulse b1 [uT]
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename

%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = defs.B1pa*gamma_rad*defs.tp; % flip angle of sat pulse
% create pulseq saturation pulse object
%satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', t_p,'system',lims,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
satPulse      = mr.makeSincPulse(fa_sat, 'Duration', defs.tp, 'system', seq.sys,'timeBwProduct', 2,'apodization', 0.15); % philips-like sinc

[B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,defs.tp,defs.td,1,gamma_hz);
defs.B1rms = B1rms;


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

%% plot
saveSaturationPhasePlot(seq_filename);

%% call standard sim
M_z = simulate_pulseqcest(seq_filename,'../../sim-library/WM_3T_default_7pool_bmsim.yaml');

%% plot
plotSimulationResults(M_z,defs.offsets_ppm, defs.M0_offset);
writematrix(M_z', ['M_z_' seq_filename '.txt']);
