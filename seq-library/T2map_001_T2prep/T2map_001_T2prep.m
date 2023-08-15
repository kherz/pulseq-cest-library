%% T2map_001_T2prep
% Creates a sequence file for T2 mapping
% T2prep is realized by three pulses of FA 90-180-90 and phases
% 90-180-(-90)
% see Figure 6 of https://doi.org/10.3348/kjr.2017.18.1.113 
% Moritz Zaiss 2023
% moritz.zaiss@fau.de

% author name for sequence file
author = 'Moritz Zaiss';

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
defs.tp            = 1.2e-3         ; % [s] % Single saturation pulse length
defs.Trec          = 10              ; % recovery time [s]
defs.TE            = [0 0.01 0.025 0.03 0.04 0.05 0.1 0.2 0.3 0.5 1.0];
defs.offsets_ppm   = zeros(1, numel(defs.TE));
defs.num_meas      = numel(defs.TE); % number of repetitions
defs.seq_id_string = seqid           ; % unique seq id
defs.B0            = -1              ; % dummy b0
defs.B1pa        = 20;  % mean sat pulse b1 [uT]
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename


%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]

rB1=1.0;  % this is just to show the effect of B1 inhomogeneities, e.g. set to rB1=0.8

% create pulseq object
pre90Pulse = mr.makeBlockPulse(rB1*pi/2, 'Duration', defs.tp , 'system',seq.sys);      % flip down (ex)
refocPulse = mr.makeBlockPulse(rB1*pi, 'Duration', defs.tp ,'Phase', pi/2,'system',seq.sys);  % refocus in direction of fan (spin-lock-like)
post90Pulse = mr.makeBlockPulse(rB1*pi/2, 'Duration', defs.tp ,'Phase', pi,'system',seq.sys); % flip up, same FA as ex, pi phase resp. ex

% spoilers
rampTime = 1e-3;
spoilDuration = 4500e-6 + rampTime; % in s

lims = seq.sys;
spoilAmplitude = 0.8 .* lims.maxGrad; % in Hz/m % FG: hard coded to 80% of maximum gradient, more should be possible
gxSpoil=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gySpoil=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gzSpoil=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);

%% loop through echo times (TE)
% init sequence


for ii = 1:numel(defs.TE)
    realTehalf = defs.TE(ii)/2-lims.rfRingdownTime-lims.rfDeadTime;
    % recover time
    
    seq.addBlock(mr.makeDelay(defs.Trec));
    if defs.TE(ii) > 0 % m0 at the beginning
        
        seq.addBlock(pre90Pulse);
        seq.addBlock(mr.makeDelay(realTehalf));
        seq.addBlock(refocPulse);
        seq.addBlock(mr.makeDelay(realTehalf));
        seq.addBlock(post90Pulse);
        
        seq.addBlock(mr.makeDelay(100e-6)); % prespoiler
        if spoiling == 1 && ii > 1 % spoil before readout after saturation
            seq.addBlock(gxSpoil,gySpoil,gzSpoil);
            seq.addBlock(mr.makeDelay(100e-6)); % necessary?
        end
    end
    
    % add pseudo adc / readout trigger event
    seq.addPseudoADCBlock();
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
figure, plot(defs.TE,M_z);
xlabel('TE [s]');
ylabel('Z [a.u.]')

writematrix(M_z', ['M_z_' seq_filename '.txt']);

% cftool(defs.TE,M_z);