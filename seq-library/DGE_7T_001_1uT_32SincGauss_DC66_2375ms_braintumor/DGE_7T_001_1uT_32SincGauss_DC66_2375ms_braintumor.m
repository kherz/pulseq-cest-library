%% DGE_7T_001_1uT_32SincGauss_DC66_2375ms_braintumor
% Creates a sequence file for a DGE protocol with Sinc-Gaussian pulses, 66% DC and tsat of 2.3 s
% Xu X, Yadav NN, Knutsson L, et al. Dynamic Glucose-Enhanced (DGE) MRI: Translation to Human Scanning and First Results in Glioma Patients. Tomography. 2015;1(2):105-114. doi:10.18383/j.tom.2015.00175
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
defs.n_pulses      = 32              ; % number of pulses
defs.tp            = 50e-3           ; % pulse duration [s]
defs.td            = 25e-3            ; % interpulse delay [s]
defs.Trec          = 2             ; % recovery time [s]
defs.DCsat         = (defs.tp)/(defs.tp+defs.td); % duty cycle
sat_offset = 1.2;
defs.num_meas      = 10 ; % should b enough for a sketch
defs.offsets_ppm   = []; % offset vector [ppm]
for nmeas = 1:defs.num_meas
    defs.offsets_ppm   = [defs.offsets_ppm sat_offset];
end

defs.Tsat          = defs.n_pulses*(defs.tp+defs.td) - ...
    defs.td ;  % saturation time [s]
defs.FREQ		   = 298.0348 ;        % Approximately 7 T 
defs.B0            = defs.FREQ/(gamma_hz);  % Calculate B0    
defs.seq_id_string = seqid           ; % unique seq id

defs.B1peak      = 1.96;  % peak b1 of saturation pulse 
defs.spoiling    = 1;     % 0=no spoiling, 1=before readout, Gradient in x,y,z

seq_filename = strcat(defs.seq_id_string,'.seq'); % filename


%% create scanner events
% satpulse
gamma_rad = gamma_hz*2*pi;        % [rad/uT]
fa_sat        = gamma_rad*defs.tp; % dummy flip angle
% create pulseq saturation pulse object

%satPulse      = mr.makeGaussPulse(fa_sat, 'Duration', t_p,'system',lims,'timeBwProduct', 0.2,'apodization', 0.5); % siemens-like gauss
satPulse      = mr.makeSincPulse(fa_sat, 'Duration', defs.tp, 'system', seq.sys,'timeBwProduct', 2,'apodization', 0.15); % philips-like sinc
satPulse.signal = (satPulse.signal./max(satPulse.signal)).*defs.B1peak*gamma_hz; 


[B1rms,B1cwae,B1cwae_pure,alpha]= calculatePowerEquivalents(satPulse,defs.tp,defs.td,1,gamma_hz);
defs.B1rms = B1rms;


%% loop through zspec offsets
offsets_Hz = defs.offsets_ppm*defs.FREQ;

% loop through offsets and set pulses and delays
for currentOffset = offsets_Hz
    
    seq.addBlock(mr.makeDelay(defs.Trec)); % recovery time
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
