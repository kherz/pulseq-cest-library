%% MRF_CEST_MT_3T_003_13SL_DC50_2500ms
% Creates a sequence file for a MT CEST Fingerprinting protocol
% Reference: Perlman, O., Ito, H., Herz, K. et al. 
% Quantitative imaging of apoptosis following oncolytic virotherapy by magnetic resonance fingerprinting aided by deep learning. 
% Nat. Biomed. Eng, 2022, https://doi.org/10.1038/s41551-021-00809-7
%

% author name for sequence file
author = 'Kai Herz';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% scanner limits
% choose less demanding rf limits for very short sl breaks
lims = getScannerLimits();
lims.rfDeadTime     = 200e-6;
lims.rfRingdownTime = 50e-6;
seq = SequenceSBB(lims);
gamma_hz  =seq.sys.gamma*1e-6;                  % for H [Hz/uT]

%% varying MRF parameters
num_meas    = 30;
TR          = ones(num_meas,1)*3.5;
Tsat        = ones(num_meas,1)*2.5;
offsets_ppm = [8,   6,   6,   10,  10,  10,  8,   6,   8,   14,  14,  10,  6,   10,  8,   6,   8,   10,  14,  14,  6,   8,   14,  6,   14,  14,  14,  8,   10,  8  ];
B1          = [2.0, 2.0, 1.7, 1.5, 1.2, 1.2, 3.0, 0.5, 3.0, 1.0, 2.2, 3.2, 1.5, 0.7, 1.5, 2.2, 2.5, 1.2, 3.0, 0.2, 1.5, 2.5, 0.7, 4.0, 3.2, 3.5, 1.5, 2.7, 0.7, 0.5];

%% sequence definitions
% everything in defs gets written as definition in .seq-file
defs.n_pulses      = 13            ; % number of pulses
defs.tp            = 100e-3        ; % pulse duration [s]
defs.td            = 100e-3        ; % interpulse delay [s]
defs.Trec          = TR - Tsat     ; % recovery time [s]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
defs.num_meas      = num_meas      ; % number of measurements
defs.offsets_ppm   = offsets_ppm   ; % offset vector [ppm]
defs.Tsat          = Tsat          ;  % saturation time [s]
defs.FREQ		   = 127.7292 ;         % Approximately 3 T 
defs.B0            = defs.FREQ/(gamma_hz);  % Calculate B0    
defs.seq_id_string = seqid         ; % unique seq id
defs.B1rms        = B1;
defs.spoiling    = 1;                    % 0=no spoiling, 1=before readout, Gradient in x,y,z
seq_filename = strcat(defs.seq_id_string,'.seq'); % filename


%% spin lock specific preparation
% td should be between block pulses, so we fit the tipping pulses in the
% delay between pulses
tp_sl        = 1e-3;                                    % duration of tipping pulse for sl
td_sl        = (lims.rfDeadTime + lims.rfRingdownTime); % delay between tip and sat pulse
slTimePerSat = 2*(tp_sl+td_sl); % additional time of sl pulses for 1 sat pulse
% check if sl pulses fit in DC
if any(defs.Trec < slTimePerSat)
    error('DC too high for SL prepatration pulses!');
end

%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
offsets_Hz = defs.offsets_ppm*defs.FREQ;

% loop through measurements
for m = 1:num_meas
    seq.addBlock(mr.makeDelay(defs.Trec(m))); % recovery time
    % calculate spin lock pulses for current B1
    cB1 = B1(m); % get current B1 from schedule
    satFa = cB1*gamma_rad*defs.tp;  % saturation pulse flip angle
    faSL = atan(cB1/(offsets_ppm(m)*defs.B0));
    
    preSL = mr.makeBlockPulse(faSL,'Duration',tp_sl, 'Phase', -pi/2,'system',lims);
    satPulse = mr.makeBlockPulse(satFa,'Duration',defs.tp, 'freqOffset', offsets_Hz(m),'system',lims);
    accumPhase = mod(offsets_Hz(m)*2*pi*defs.tp,2*pi);
    postSL = mr.makeBlockPulse(faSL,'Duration',tp_sl, 'Phase', accumPhase+pi/2,'system',lims);
    satPulse.freqOffset = offsets_Hz(m); % set freuqncy offset of the pulse
    for np = 1:defs.n_pulses
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
        if np < defs.n_pulses % delay between pulses
            seq.addBlock(mr.makeDelay(defs.td-slTimePerSat)); % add delay
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