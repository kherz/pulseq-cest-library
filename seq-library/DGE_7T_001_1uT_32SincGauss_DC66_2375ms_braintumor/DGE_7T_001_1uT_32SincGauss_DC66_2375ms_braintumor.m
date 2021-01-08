%% DGE_7T_001_1uT_32SincGauss_DC66_2375ms_braintumor
% Creates a sequence file for a DGE protocol with Sinc-Gaussian pulses, 66% DC and tsat of 2.3 s
% Xu X, Yadav NN, Knutsson L, et al. Dynamic Glucose-Enhanced (DGE) MRI: Translation to Human Scanning and First Results in Glioma Patients. Tomography. 2015;1(2):105-114. doi:10.18383/j.tom.2015.00175
%
% Kai Herz 2020
% kai.herz@tuebingen.mpg.de

% author name for sequence file
author = 'Kai Herz';

%% get id of generation file
if strcmp(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    [~, seqid] = fileparts(matlab.desktop.editor.getActiveFilename);
else
    [~, seqid] = fileparts(which(mfilename));
end

%% sequence definitions
% everything in seq_defs gets written as definition in .seq-file
seq_defs.n_pulses      = 32              ; % number of pulses
seq_defs.tp            = 50e-3           ; % pulse duration [s]
seq_defs.td            = 25e-3            ; % interpulse delay [s]
seq_defs.Trec          = 2             ; % recovery time [s]
seq_defs.DCsat         = (seq_defs.tp)/(seq_defs.tp+seq_defs.td); % duty cycle
sat_offset = 1.2;
seq_defs.num_meas      = 10 ; % should b enough for a sketch
seq_defs.offsets_ppm   = []; % offset vector [ppm]
for nmeas = 1:seq_defs.num_meas
    seq_defs.offsets_ppm   = [seq_defs.offsets_ppm sat_offset];
end

seq_defs.Tsat          = seq_defs.n_pulses*(seq_defs.tp+seq_defs.td) - ...
    seq_defs.td ;  % saturation time [s]
seq_defs.B0            = 7               ; % B0 [T]
seq_defs.seq_id_string = seqid           ; % unique seq id


%% get info from struct
offsets_ppm = seq_defs.offsets_ppm; % [ppm]
Trec        = seq_defs.Trec;        % recovery time between scans [s]

tp          = seq_defs.tp;          % sat pulse duration [s]
td          = seq_defs.td;          % delay between pulses [s]
n_pulses    = seq_defs.n_pulses;    % number of sat pulses per measurement. if DC changes use: n_pulses = round(2/(t_p+t_d))
B0          = seq_defs.B0;          % B0 [T]
B1peak      = 1.96;  % peak b1 of saturation pulse 
spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(seq_defs.seq_id_string,'.seq'); % filename

%% scanner limits
% see pulseq doc for more ino
lims = Get_scanner_limits();

%% create scanner events
% satpulse
gyroRatio_hz  = 42.5764;                  % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;        % [rad/uT]
fa_sat        = gyroRatio_rad*tp; % dummy flip angle
% create pulseq saturation pulse object

%satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', t_p,'system',lims,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
satPulse      = mr.makeSincPulse(fa_sat, 'Duration', tp, 'system', lims,'timeBwProduct', 2,'apodization', 0.15); % philips-like sinc
satPulse.signal = (satPulse.signal./max(satPulse.signal)).*B1peak*gyroRatio_hz; 


[B1cwpe,B1cwae,B1cwae_pure,alpha]= calc_power_equivalents(satPulse,tp,td,1,gyroRatio_hz);
seq_defs.B1cwpe = B1cwpe;

% spoilers
spoilRiseTime = 1e-3;
spoilDuration = 4500e-6+ spoilRiseTime; % [s]
% create pulseq gradient object
[gxSpoil, gySpoil, gzSpoil] = Create_spoiler_gradients(lims, spoilDuration, spoilRiseTime);

% pseudo adc, not played out
pseudoADC = mr.makeAdc(1,'Duration', 1e-3);

%% loop through zspec offsets
offsets_Hz = offsets_ppm*gyroRatio_hz*B0;

% init sequence
seq = mr.Sequence();
% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    
    seq.addBlock(mr.makeDelay(Trec)); % recovery time
    satPulse.freqOffset = currentOffset; % set freuqncy offset of the pulse
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
    if spoiling % spoiling before readout
        seq.addBlock(gxSpoil,gySpoil,gzSpoil);
    end
    seq.addBlock(pseudoADC); % readout trigger event
end



%% write definitions
def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_filename, author);

%% plot
Run_pulseq_cest_Simulation(seq_filename,'GM_3T_001_bmsim.yaml',1,1);




