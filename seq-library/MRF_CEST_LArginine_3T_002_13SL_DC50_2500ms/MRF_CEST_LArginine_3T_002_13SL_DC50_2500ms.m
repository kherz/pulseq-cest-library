%% MRF_CEST_L-Arginine_3T_002_13SL_DC50_2500ms
% Creates a sequence file for a L-Arginine CEST Fingerprinting protocol
% Reference: tbd
%

% author name for sequence file
author = 'Kai Herz';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% varying MRF parameters
num_meas    = 30;
TR          = ones(num_meas,1)*3.5;
Tsat        = ones(num_meas,1)*2.5;
offsets_ppm = ones(num_meas,1)*3.0;
B1          = [2 2 1.7 1.5 1.2 1.2 3 0.5 3 1 2.2 3.2 1.5 0.7 1.5 2.2 2.5 1.2 3 0.2 1.5 2.5 0.7 4 3.2 3.5 1.5 2.7 0.7 0.5];

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = 13            ; % number of pulses
seq_defs.tp            = 100e-3        ; % pulse duration [s]
seq_defs.td            = 100e-3        ; % interpulse delay [s]
seq_defs.Trec          = TR - Tsat     ; % recovery time [s]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.num_meas      = num_meas      ; % number of measurements
seq_defs.offsets_ppm   = offsets_ppm   ; % offset vector [ppm]
seq_defs.Tsat          = Tsat          ;  % saturation time [s]
seq_defs.B0            = 3             ; % B0 [T]
seq_defs.seq_id_string = seqid         ; % unique seq id
seq_defs.B1rms        = B1;


%% get info from struct
Trec        = seq_defs.Trec;        % recovery time between scans [s]
tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B0          = seq_defs.B0;          % B0 [T]
spoiling    = 1;                    % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% choose less demanding rf limits for very short sl breaks
lims = getScannerLimits();
lims.rfDeadTime     = 200e-6;
lims.rfRingdownTime = 50e-6;
seq = SequenceSBB(lims);

%% spin lock specific preparation
% td should be between block pulses, so we fit the tipping pulses in the
% delay between pulses
tp_sl        = 1e-3;                                    % duration of tipping pulse for sl
td_sl        = (lims.rfDeadTime + lims.rfRingdownTime); % delay between tip and sat pulse
slTimePerSat = 2*(tp_sl+td_sl); % additional time of sl pulses for 1 sat pulse
% check if sl pulses fit in DC
if any(Trec < slTimePerSat)
    error('DC too high for SL prepatration pulses!');
end

%% create scanner events
% satpulse
gyroRatio_hz  = 42.5764;                  % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;        % [rad/uT]
offsets_Hz = offsets_ppm*gyroRatio_hz*B0;

% loop through measurements
for m = 1:num_meas
    seq.addBlock(mr.makeDelay(Trec(m))); % recovery time
    % calculate spin lock pulses for current B1
    cB1 = B1(m); % get current B1 from schedule
    satFa = cB1*gyroRatio_rad*tp;  % saturation pulse flip angle
    faSL = atan(cB1/(offsets_ppm(m)*B0));
    
    preSL = mr.makeBlockPulse(faSL,'Duration',tp_sl, 'Phase', -pi/2,'system',lims);
    satPulse = mr.makeBlockPulse(satFa,'Duration',tp, 'freqOffset', offsets_Hz(m),'system',lims);
    accumPhase = mod(offsets_Hz(m)*2*pi*tp,2*pi);
    postSL = mr.makeBlockPulse(faSL,'Duration',tp_sl, 'Phase', accumPhase+pi/2,'system',lims);
    satPulse.freqOffset = offsets_Hz(m); % set freuqncy offset of the pulse
    for np = 1:n_pulses
        % add sl phase cycling. We dont need accumulated phse here, as z
        % mag gets restored after every sat pulse
        phaseCycling         =  50/180*pi;
        preSL.phaseOffset    = preSL.phaseOffset    + phaseCycling;
        postSL.phaseOffset   = postSL.phaseOffset   + phaseCycling;
        satPulse.phaseOffset = satPulse.phaseOffset + phaseCycling;
        % add pulses
        seq.addBlock(preSL);
        seq.addBlock(satPulse); 
        seq.addBlock(postSL);
        % calc phase for next rf pulse
        if np < n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(td-slTimePerSat)); % add delay
        end
    end
    if spoiling % spoiling before readout
        seq.addSpoilerGradients()
    end
    seq.addPseudoADCBlock(); % readout trigger event
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