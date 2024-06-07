%% WASABIT1_3T_001_3p7uT_1block_5p12ms
% Creates a sequence file for a WASABITI protocol
%

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
defs.n_pulses      = 1; % number of pulses
defs.tp            = 5.12e-3         ; % pulse duration [s]
defs.Trec          = [15 0.5 1 2 3 4 5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 5 4 3 2 1 0.5]              ; % recovery time [s]
defs.offsets_ppm   = [-300 -4 -3.73333333 -3.46666667 -3.2 -2.93333333 -2.66666667 -2.4 -2.13333333 -1.86666667 -1.6 -1.33333333 -1.06666667 -0.8 -0.533333333 -0.266666667 0 0.266666667 0.533333333 0.8 1.06666667 1.33333333 1.6 1.86666667 2.13333333 2.4 2.66666667 2.93333333 3.2 3.46666667 3.73333333 4 ];
defs.num_meas      = numel(defs.offsets_ppm); % number of repetition
defs.seq_id_string = seqid           ; % unique seq id
defs.FREQ		   = 127.7292;          % Approximately 3 T 
defs.B1pa        = 3.7;  % mean sat pulse b1 [uT]
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename


%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]

fa_sat        = defs.B1pa*gamma_rad*defs.tp; % flip angle of sat pulse

% create pulseq saturation pulse object
satPulse      = mr.makeBlockPulse(fa_sat, 'Duration', defs.tp, 'system', seq.sys); % block pulse


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


%% loop through rec times and offsets
% init sequence
hs_pulse.freqOffset = 0;

ii=1;
for t_rec = defs.Trec
    
    for np = 1:3
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
    seq.addBlock(mr.makeDelay(t_rec));
    
     % add single pulse
    satPulse.freqOffset = defs.offsets_ppm(ii)*defs.FREQ; % set freuqncy offset of the pulse
    seq.addBlock(satPulse) % add sat pulse
    if defs.spoiling % spoiling before readout
        seq.addSpoilerGradients();
    end
    ii=ii+1;
    
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
figure, plot(M_z);
xlabel('TI [s]');
ylabel('Z [a.u.]')


writematrix(M_z', ['M_z_' seq_filename '.txt']);