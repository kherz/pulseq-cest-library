%% T1map_001_3HypSec
% Creates a sequence file for a T1map protocol
%
% Patrick Schuenke 2021
% patrick.schuenke@ptb.de

% author name for sequence file
author = 'Patrick Schuenke';

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
defs.tp            = 8e-3         ; % pulse duration [s]
defs.Trec          = 1              ; % recovery time [s]
defs.TI            = [10 6 5 4 3 2.5 2 1.5 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1];
defs.offsets_ppm   = zeros(1, numel(defs.TI));
defs.num_meas      = numel(defs.offsets_ppm); % number of repetition
defs.seq_id_string = seqid           ; % unique seq id
defs.B0            = -1              ; % dummy b0
defs.B1pa        = 20;  % mean sat pulse b1 [uT]
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename


%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]

% create pulseq object
hs_pulse = makeHSHalfPassagePulse(20,seq.sys);
defs.B1rms = defs.B1pa;

% spoilers
rampTime = 1e-3;
spoilDuration = 4500e-6 + rampTime; % in s

lims = seq.sys;
spoilAmplitude = 0.8 .* lims.maxGrad; % in Hz/m % FG: hard coded to 80% of maximum gradient, more should be possible
gxSpoil1=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gySpoil1=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gzSpoil1=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);

spoilAmplitude = -0.7 .* lims.maxGrad;
gxSpoil2=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gySpoil2=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gzSpoil2=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);

spoilAmplitude = 0.6 .* lims.maxGrad;
gxSpoil3=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gySpoil3=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);
gzSpoil3=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration, 'riseTime', rampTime, 'system',lims);


%% loop through "inversion" times (TI)
% init sequence
hs_pulse.freqOffset = 0;

for t_prep = defs.TI
    if defs.Trec > 0
        seq.addBlock(mr.makeDelay(defs.Trec)); % recovery time
    end
        
    for np = 1:defs.n_pulses
        seq.addBlock(hs_pulse);
        if mod(np-1,3) == 0
            seq.addBlock(gxSpoil1,gySpoil2,gzSpoil3);
        elseif mod(np-1,2) == 0
            seq.addBlock(gxSpoil3,gySpoil1,gzSpoil2);
        else
            seq.addBlock(gxSpoil2,gySpoil3,gzSpoil1);
        end
    end
    
    % add variable "inversion" time delay
    seq.addBlock(mr.makeDelay(t_prep));
    
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
figure, plot(defs.TI,M_z);
xlabel('TI [s]');
ylabel('Z [a.u.]')


writematrix(M_z', ['M_z_' seq_filename '.txt']);