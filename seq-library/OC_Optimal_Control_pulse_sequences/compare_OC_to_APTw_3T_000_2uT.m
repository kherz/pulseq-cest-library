%% Compare OC pulses to APTw_3T_000_2uT_1block_2s_braintumor standard
% Generates OC pulse train with matched B1rms and saturation time
% Clemens Stilianu 2025, stilianu@hotmail.com

% This script compares OC pulse performance against the APTw brain tumor
% standard (APTw_3T_000_2uT_1block_2s_braintumor). The pulse train uses
% 18x100ms pulses with 90% duty cycle to match the 2s CW saturation time.
% Simulation uses WM_3T_default_7pool_bmsim.yaml for direct comparison
% with the reference spectrum M_z_APTw_3T_000_2uT_1block_2s_braintumor.seq.txt.

%% Parameters matched to standard
seq = SequenceSBB(getScannerLimits());
gamma_hz = seq.sys.gamma * 1e-6;

% Pulse train parameters for exactly 2s total with 90% DC
n_pulses = 18;              % 18 pulses
tp = 100e-3;                % 100ms each = 1.8s pulse time
Tsat_target = 2.0;          % target total saturation time [s]

% Calculate td to get exactly 2s total (train ends with pulse, not delay)
td = (Tsat_target - n_pulses * tp) / (n_pulses - 1);  % ~11.76ms

% Calculate duty cycle and B1 scaling
Tsat = n_pulses * tp + (n_pulses - 1) * td;     % total saturation time = 2s
duty_cycle = (n_pulses * tp) / Tsat;            % 90%

% Scale B1 so train B1rms = 2uT (matching standard)
B1_train_target = 2.0;                          % uT (same as standard)
B1_pulse = B1_train_target / sqrt(duty_cycle);  % ~2.11 uT

% Sequence definitions
defs.n_pulses = n_pulses;
defs.tp = tp;
defs.td = td;
defs.Trec = 3.5;
defs.Trec_M0 = 3.5;
defs.M0_offset = -300;
defs.offsets_ppm = [defs.M0_offset -4:0.25:4];  % match standard
defs.FREQ = 127.7292;                            % match standard 3T
defs.B0 = defs.FREQ / gamma_hz;
defs.B1pa = B1_pulse;
defs.spoiling = 1;
defs.DCsat = duty_cycle;
defs.Tsat = Tsat;
defs.num_meas = numel(defs.offsets_ppm);
defs.seq_id_string = 'OC_APTw_3T_000_2uT_18pulses';

seq_filename = [defs.seq_id_string '.seq'];

%% Generate OC pulse
fa_sat = B1_pulse * gamma_hz * 2 * pi * tp;
satPulse = makeOCPulse(fa_sat, 'system', seq.sys, 'duration', tp, 'useLowPass', true);

% Get initial B1rms (calculatePowerEquivalents returns train RMS including DC)
[B1rms_init, ~, ~, ~] = calculatePowerEquivalents(satPulse, tp, td, 1, gamma_hz);

% Scale pulse so train B1rms = 2 uT
scale_factor = B1_train_target / B1rms_init;
satPulse.signal = satPulse.signal * scale_factor;

% Verify final power equivalents
[B1rms, B1cwae, ~, ~] = calculatePowerEquivalents(satPulse, tp, td, 1, gamma_hz);
B1rms_pulse = rms(satPulse.signal / gamma_hz);
defs.B1rms = B1rms;

fprintf('OC pulse train characteristics:\n');
fprintf('  Pulses: %d x %.0f ms\n', n_pulses, tp*1000);
fprintf('  Interpulse delay: %.2f ms\n', td*1000);
fprintf('  Total duration: %.0f ms\n', Tsat*1000);
fprintf('  Duty cycle: %.1f%%\n', duty_cycle*100);
fprintf('  B1 pulse RMS: %.2f uT\n', B1rms_pulse);
fprintf('  B1 train RMS: %.2f uT\n', B1rms);
fprintf('  B1 CWAE: %.2f uT\n', B1cwae);

%% Build sequence
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

    for np = 1:n_pulses
        satPulse.phaseOffset = mod(accumPhase, 2*pi);
        seq.addBlock(satPulse);
        accumPhase = mod(accumPhase + currentOffset*2*pi*(numel(find(abs(satPulse.signal)>0))*1e-6), 2*pi);
        if np < n_pulses
            seq.addBlock(mr.makeDelay(td));
        end
    end

    if defs.spoiling
        seq.addSpoilerGradients();
    end
    seq.addPseudoADCBlock();
end

%% Write sequence
def_fields = fieldnames(defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, defs.(def_fields{n_id}));
end
seq.write(seq_filename, 'Clemens Stilianu');
fprintf('Sequence written: %s\n', seq_filename);

%% Simulate OC sequence
yaml_file = '../../sim-library/WM_3T_default_7pool_bmsim.yaml';
M_z_oc = simulate_pulseqcest(seq_filename, yaml_file);

%% Load standard result
standard_file = '../APTw_3T_000_2uT_1block_2s_braintumor/M_z_APTw_3T_000_2uT_1block_2s_braintumor.seq.txt';
M_z_std = readmatrix(standard_file);

%% Plot comparison
offsets = defs.offsets_ppm;
offsets_plot = offsets(2:end);  % exclude M0

figure('Name', 'OC vs Standard APTw Comparison');

% Z-spectrum
subplot(1, 2, 1);
plot(offsets_plot, M_z_std(2:end), 'k-', 'LineWidth', 2); hold on;
plot(offsets_plot, M_z_oc(2:end), 'b--', 'LineWidth', 2);
set(gca, 'XDir', 'reverse');
xlabel('Offset [ppm]');
ylabel('M_z/M_0');
legend({'Standard (2s block)', 'OC (18x100ms)'}, 'Location', 'southeast');
title('Z-spectrum');
grid on;

% MTRasym
subplot(1, 2, 2);
center_idx = find(offsets_plot == 0);
if isempty(center_idx)
    [~, center_idx] = min(abs(offsets_plot));
end
pos_idx = center_idx+1:length(offsets_plot);

asym_std = zeros(size(pos_idx));
asym_oc = zeros(size(pos_idx));

for i = 1:length(pos_idx)
    pos_offset = offsets_plot(pos_idx(i));
    neg_idx = find(abs(offsets_plot + pos_offset) < 0.01, 1);
    if ~isempty(neg_idx)
        asym_std(i) = M_z_std(neg_idx+1) - M_z_std(pos_idx(i)+1);
        asym_oc(i) = M_z_oc(neg_idx+1) - M_z_oc(pos_idx(i)+1);
    end
end

plot(offsets_plot(pos_idx), asym_std*100, 'k-', 'LineWidth', 2); hold on;
plot(offsets_plot(pos_idx), asym_oc*100, 'b--', 'LineWidth', 2);
set(gca, 'XDir', 'reverse');
xlabel('Offset [ppm]');
ylabel('MTR_{asym} [%]');
legend({'Standard (2s block)', 'OC (18x100ms)'}, 'Location', 'northeast');
title('MTR asymmetry');
grid on;

% Save result
writematrix(M_z_oc', ['M_z_' seq_filename '.txt']);
fprintf('Results saved: M_z_%s.txt\n', seq_filename);
