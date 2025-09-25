%% Oprimal Control Pulse seq generation and comparison
% Generate and compare OC vs Gaussian pulses with matched B1 RMS for pulse
% durations >= 50 ms
% Clemens Stilianu 2024, stilianu@hotmail.com

%% Parameters
B1_target = 1.0;              % µT RMS
tp = 100e-3;                  % pulse duration [s]
n_pulses = 10;                 % pulse train length
td = 10e-3;                   % interpulse delay [s]
offsets_ppm = -5:0.1:5;      % CEST offsets
sim_file = 'OC_APT_phantom_sim_test.yaml'; % simulation yaml

seq = SequenceSBB(getScannerLimits());
gamma_hz = seq.sys.gamma * 1e-6;
fa_sat = B1_target * gamma_hz * 2 * pi * tp;

defs.n_pulses = n_pulses;
defs.tp = tp;
defs.td = td;
defs.Trec = 3.5;
defs.Trec_M0 = 3.5;
defs.M0_offset = -300;
defs.offsets_ppm = offsets_ppm;
defs.FREQ = 123.2627;
defs.B0 = defs.FREQ/(gamma_hz);
defs.B1pa = B1_target;
defs.spoiling = 1;

%% Generate OC sequence
satPulse_oc = makeOCPulse(fa_sat, 'system', seq.sys, 'duration', tp, 'useLowPass', true);

B1_oc = satPulse_oc.signal / gamma_hz;
RMS_oc = rms(B1_oc);
satPulse_oc.signal = satPulse_oc.signal * (B1_target / RMS_oc);

defs.seq_id_string = 'OC_minimal';
generateSequence(satPulse_oc, defs);
saveSaturationPhasePlot('OC_minimal.seq');

%% Generate Gaussian sequence  
satPulse_gauss = mr.makeGaussPulse(fa_sat, 'Duration', tp, 'system', seq.sys, ...
                                  'timeBwProduct', 0.2, 'apodization', 0.5);

B1_gauss = satPulse_gauss.signal / gamma_hz;
RMS_gauss = rms(B1_gauss);
satPulse_gauss.signal = satPulse_gauss.signal * (B1_target / RMS_gauss);

defs.seq_id_string = 'Gauss_minimal';
generateSequence(satPulse_gauss, defs);
saveSaturationPhasePlot('Gauss_minimal.seq');

%% Plot pulses
t_oc = (0:length(B1_oc)-1) * seq.sys.rfRasterTime * 1000;
t_gauss = (0:length(B1_gauss)-1) * seq.sys.rfRasterTime * 1000;

B1_oc_final = satPulse_oc.signal / gamma_hz;
B1_gauss_final = satPulse_gauss.signal / gamma_hz;

figure;
plot(t_oc, B1_oc_final, 'b-', 'LineWidth', 2); hold on;
plot(t_gauss, B1_gauss_final, 'r-', 'LineWidth', 2);
xlabel('Time [ms]'); ylabel('B1 [µT]'); legend({'OC', 'Gaussian'});

%% Simulate and plot results
if exist(sim_file, 'file')
    M_oc = simulate_pulseqcest('OC_minimal.seq', sim_file);
    M_gauss = simulate_pulseqcest('Gauss_minimal.seq', sim_file);
    
    % MTR asymmetry
    center_idx = find(offsets_ppm == 0);
    if isempty(center_idx), [~, center_idx] = min(abs(offsets_ppm)); end
    pos_idx = center_idx+1:length(offsets_ppm);
    
    asym_oc = zeros(size(pos_idx));
    asym_gauss = zeros(size(pos_idx));
    
    for i = 1:length(pos_idx)
        pos_offset = offsets_ppm(pos_idx(i));
        neg_idx = find(abs(offsets_ppm + pos_offset) < 0.01, 1);
        if ~isempty(neg_idx)
            asym_oc(i) = M_oc(neg_idx) - M_oc(pos_idx(i));
            asym_gauss(i) = M_gauss(neg_idx) - M_gauss(pos_idx(i));
        end
    end
    
    % Combined CEST spectra plot
    figure;
    subplot(1,2,1);
    plot(offsets_ppm, M_oc, 'bo-', 'LineWidth', 2); hold on;
    plot(offsets_ppm, M_gauss, 'rs-', 'LineWidth', 2);
    xlabel('Offset [ppm]'); ylabel('Mz/M0'); legend({'OC', 'Gaussian'});
    set(gca, 'XDir', 'reverse');
    
    subplot(1,2,2);
    plot(offsets_ppm(pos_idx), asym_oc*100, 'bo-', 'LineWidth', 2); hold on;
    plot(offsets_ppm(pos_idx), asym_gauss*100, 'rs-', 'LineWidth', 2);
    xlabel('Offset [ppm]'); ylabel('MTR asym [%]'); legend({'OC', 'Gaussian'});
    set(gca, 'XDir', 'reverse');
else
    fprintf('No simulation file found: %s\n', sim_file);
end

%% Helper function
function generateSequence(satPulse, defs)
    seq = SequenceSBB(getScannerLimits());
    offsets_Hz = defs.offsets_ppm * defs.FREQ;
    
    for currentOffset = offsets_Hz
        if currentOffset == defs.M0_offset * defs.FREQ
            if defs.Trec_M0 > 0
                seq.addBlock(mr.makeDelay(defs.Trec_M0));
            end
        else
            if defs.Trec > 0
                seq.addBlock(mr.makeDelay(defs.Trec));
            end
        end
        
        satPulse.freqOffset = currentOffset;
        accumPhase = 0;
        
        for np = 1:defs.n_pulses
            satPulse.phaseOffset = mod(accumPhase, 2*pi);
            seq.addBlock(satPulse);
            accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6), 2*pi);
            if np < defs.n_pulses
                seq.addBlock(mr.makeDelay(defs.td));
            end
        end
        
        if defs.spoiling
            seq.addSpoilerGradients();
        end
        seq.addPseudoADCBlock();
    end
    
    def_fields = fieldnames(defs);
    for n_id = 1:numel(def_fields)
        seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
    end
    seq.write([defs.seq_id_string '.seq'], 'Clemens Stilianu');
end