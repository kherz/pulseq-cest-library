%% APTw_3T_004_2uT_19SL_DC95_2s_braintumor
% this file requires pulseq-cest sim commit
% 94484e494e4897287ab0c78160bf654106feeff7 or later, as previously weff sign  error lead to artifacts
% see also https://github.com/kherz/pulseq-cest/pull/6
%
%
% Creates a sequence file for an APTw protocol with conventional offresonant SL pulses, 50% DC and tsat of 2 s
% This is a pulseq spin-lock pulse train following the paper of Roellofs et al.
% citation:
% Roeloffs, V., Meyer, C., Bachert, P., and Zaiss, M. (2014) 
% Towards quantification of pulsed spinlock and CEST at clinical MR scanners: 
% an analytical interleaved saturation–relaxation (ISAR) approach, 
% NMR Biomed., 28, 40– 53, doi: 10.1002/nbm.3192. 

% Moritz Zaiss 2021
% moritz.zaiss@tuebingen.mpg.de

% author name for sequence file
author = 'Moritz Zaiss';

%% get id of generation file
if contains(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = [19]       ; % number of pulses
seq_defs.tp            = 100e-3          ; % pulse duration [s]
seq_defs.td            = 5e-3           ; % interpulse delay [s]
seq_defs.Trec          = 3.5             ; % approx [s]
seq_defs.M0_offset     = -300            ; % m0 offset [ppm]
seq_defs.DCsat         = seq_defs.tp/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = [seq_defs.M0_offset -4:0.1:4]; % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ; % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses.*(seq_defs.tp+ seq_defs.td)-seq_defs.td;
seq_defs.FREQ		   = 127.7292          % Approximately 3 T  
seq_defs.seq_id_string = seqid           ; % unique seq id
seq_defs.B1pa           = [2]   ;


%% get info from struct
m0_offset=seq_defs.M0_offset;
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
Trec        = seq_defs.Trec;        % recovery time between scans [s]
Trec_M0     = seq_defs.Trec;        % recovery time between scans [s]
tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B1          = seq_defs.B1pa;          % mean sat pulse b1 [uT]
spoiling    = 1;                   % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
seq = SequenceSBB(getScannerLimits());

%% create scanner events
% satpulse
gamma_hz  =seq.sys.gamma*10e-6;                  % for H [Hz/uT]
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat = B1*gamma_rad*tp;  % saturation pulse flip angle
satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', tp, 'system', seq.sys); % block pusle cw
[B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,tp,td,1,gamma_hz);

seq_defs.B1rms = B1rms;


DC=tp/(tp+td);%Duty cycle

% SL pulses are played out during the interpulse delay. We need to make
% sure that they fit in the DC here
slPrepPulseTime = 1e-3; %Gaps between pulses [s]
slPauseTime = 250e-6;


additionalSLPrepTime = 2*slPrepPulseTime+2*slPauseTime;
td = td - additionalSLPrepTime; % DC is between block bulses SL is in between
if td < 100e-6
    error('DC too high for SL prepatration pulses!');
end

% these rf times getting added in the run function of pulseq! If we want
% the timing to be exact we have to take care of this
if(slPauseTime-seq.sys.rfDeadTime-seq.sys.rfRingdownTime) <= 0
    error('slPauseTime is too short for hardware limits');
end
slPauseTime = slPauseTime-seq.sys.rfDeadTime-seq.sys.rfRingdownTime;


minFa = 0.38; % this is the limit for prep pulses to be played out


%% loop through zspec offsets
offsets_Hz = offsets_ppm*seq_defs.FREQ;


% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    if currentOffset == seq_defs.M0_offset*seq_defs.FREQ
        if Trec_M0 > 0
            seq.addBlock(mr.makeDelay(Trec_M0));
        end
    else
        if Trec > 0
            seq.addBlock(mr.makeDelay(Trec)); % recovery time
        end
    end
    
    fa_sat        = B1*gamma_rad*tp; % flip angle of sat pulse
    faSL = atan(gamma_hz*B1/(currentOffset));   % thats the angle theta of the effective system
    preSL = mr.makeBlockPulse(faSL,'Duration',slPrepPulseTime, 'Phase', -pi/2,'system',seq.sys);
    satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', tp,'freqOffset', currentOffset, 'system', seq.sys);
    accumPhase = currentOffset*360*tp*pi/180; % scanner needs the correct phase at the end of the sturation pulse
    postSL = mr.makeBlockPulse(faSL,'Duration',slPrepPulseTime, 'Phase', accumPhase+pi/2,'system',seq.sys);
        
    for np = 1:n_pulses
              
        if fa_sat == 0 % pulses with amplitude/FA 0 have to be replaced by delay
            seq.addBlock(mr.makeDelay(tp));
            seq.addBlock(mr.makeDelay(additionalSLPrepTime));
        else
            if abs(faSL) > minFa/180*pi
                seq.addBlock(preSL);
                seq.addBlock(mr.makeDelay(slPauseTime));
            else
                 seq.addBlock(mr.makeDelay(mr.calcDuration(preSL)+slPauseTime));
            end
            seq.addBlock(satPulse);
            if abs(faSL) > minFa/180*pi
                seq.addBlock(mr.makeDelay(slPauseTime));
                seq.addBlock(postSL);
            else
                seq.addBlock(mr.makeDelay(mr.calcDuration(postSL)+slPauseTime));
            end
        end
        
        if np < n_pulses % delay between pulses
            
            seq.addBlock(mr.makeDelay(td)); % add delay
            
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

%% plot
plotSimulationResults(M_z,offsets_ppm, seq_defs.M0_offset);
