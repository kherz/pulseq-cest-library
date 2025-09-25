% Clemens Stilianu 2024
% stilianu@tugraz.at 
% Implementation of generalized pulses desgined by Optimal Control (OC) pulses for pulseq CEST
% The pulse uses an arbitrary waveform loaded from a text file, providing flexibility for new experimental protocols.
% The module is intended to replace standard pulses like 'makeSincPulse' while maintaining compatibility with Pulseq sequences.
% Additonal options: 
% time interpolation: 50-75ms uses 50ms base, >=75ms uses 100ms base

function [rf, gz, gzr] = makeOCPulse(flip,varargin)
%   rf=makeOCPulse(flip, 'Duration', dur) Create OC pulse
%   with given flip angle (rad) and a fixed or user-specified duration.
%
%   rf=makeOCPulse(..., 'FreqOffset', f,'PhaseOffset',p)
%   Create OC pulse with frequency offset (Hz) and phase offset (rad).
%
%   [rf, gz]=makeOCPulse(...,'SliceThickness',st) Return the
%   slice select gradient corresponding to given slice thickness (m).
%
%   [rf, gz, gzr]=makeOCPulse(flip,lims,...) Create slice selection and 
%   slice refocusing gradients with the specificed gradient limits 
%   (e.g. amplitude, slew) and taking into account 'centerpos' parameter
%
%   See also  Sequence.addBlock

validPulseUses = {'excitation','refocusing','inversion'};

persistent parser
if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'makeOCPulse';
    
    % RF params
    addRequired(parser, 'flipAngle', @isnumeric);
    addOptional(parser, 'system', mr.opts(), @isstruct);
    addParamValue(parser, 'freqOffset', 0, @isnumeric);
    addParamValue(parser, 'phaseOffset', 0, @isnumeric);
    addParamValue(parser, 'centerpos', 0.5, @isnumeric);
    addParamValue(parser, 'duration', 100e-3, @isnumeric); % Default duration is 100 ms
    addParamValue(parser, 'useLowPass', false, @islogical); % Use low-pass filtered pulse
    % Slice params
    addParamValue(parser, 'maxGrad', 0, @isnumeric);
    addParamValue(parser, 'maxSlew', 0, @isnumeric);
    addParamValue(parser, 'sliceThickness', 0, @isnumeric);
    addParamValue(parser, 'delay', 0, @isnumeric);
    % whether it is a refocusing pulse (for k-space calculation)
    addOptional(parser, 'use', '', @(x) any(validatestring(x,validPulseUses)));
end
parse(parser, flip, varargin{:});
opt = parser.Results;

% Load the appropriate waveform based on user input
if opt.duration < 50e-3
    error('Duration must be >= 50 ms');
elseif opt.duration <= 75e-3
    % Use 50ms pulse as base (will be stretched if duration > 50ms)
    if opt.useLowPass
        pulse_data = load('OC_generalized_367_50ms_filtered.mat');
    else
        pulse_data = load('OC_generalized_367_50ms.mat');
    end
else
    % Use 100ms pulse as base (will be stretched or shrunk)
    if opt.useLowPass
        pulse_data = load('OC_generalized_370_filtered.mat');
    else
        pulse_data = load('OC_generalized_370.mat');
    end
end

field_names = fieldnames(pulse_data);
pulse_data = pulse_data.(field_names{1}); % Assuming the pulse is stored in the first field
pulse_data = pulse_data(:); % Ensure it is a column vector

% Original time vector for the loaded pulse (1e-4 s raster time)
t_original = (0:length(pulse_data)-1) * 1e-4;
original_duration = t_original(end) + 1e-4; % Total duration of original pulse

% Duration and scaling
requested_duration = opt.duration;

% Target time vector for interpolation (Pulseq raster time, 1e-6 s)
rf_raster_time = opt.system.rfRasterTime; % Typically 1e-6 s
t_target = (0:rf_raster_time:requested_duration-rf_raster_time);

% Scale the original time vector to match requested duration for proper stretching
t_original_scaled = t_original * (requested_duration / original_duration);
        
        % Debugging code removed

% Perform interpolation to match Pulseq raster time with proper stretching
pulse_data_interpolated = interp1(t_original_scaled, pulse_data, t_target, 'nearest', 'extrap'); % Use nearest, otherwise scanner wont play it
        
        % Handle any NaN or Inf values
        if any(isnan(pulse_data_interpolated)) || any(isinf(pulse_data_interpolated))
            pulse_data_interpolated(isnan(pulse_data_interpolated) | isinf(pulse_data_interpolated)) = 0;
        end

% Scale the waveform to achieve the desired flip angle
        
        % Debugging plots removed

flip_integral = sum(pulse_data_interpolated) * rf_raster_time;
if flip_integral == 0
    error('Flip integral is zero, check pulse data or interpolation.');
end
signal = pulse_data_interpolated * opt.flipAngle / (flip_integral * 2 * pi);

        

rf.type = 'rf';
rf.signal = signal;
rf.t = t_target;
rf.freqOffset = opt.freqOffset;
rf.phaseOffset = opt.phaseOffset;
rf.deadTime = opt.system.rfDeadTime;
rf.ringdownTime = opt.system.rfRingdownTime;
rf.delay = opt.delay;
if ~isempty(opt.use)
    rf.use = opt.use;
end
if rf.deadTime > rf.delay
    rf.delay = rf.deadTime;
end

if nargout > 1
    assert(opt.sliceThickness > 0,'SliceThickness must be provided');
    if opt.maxGrad > 0
        opt.system.maxGrad = opt.maxGrad;
    end
    if opt.maxSlew > 0
        opt.system.maxSlew = opt.maxSlew;
    end
    
    BW = 1 / (requested_duration / length(pulse_data_interpolated)); % Bandwidth estimation based on waveform length and requested duration
    amplitude = BW / opt.sliceThickness;
    area = amplitude * requested_duration;
    gz = mr.makeTrapezoid('z', opt.system, 'flatTime', requested_duration, 'flatArea', area);
    gzr = mr.makeTrapezoid('z', opt.system, 'Area', -area * (1 - opt.centerpos) - 0.5 * (gz.area - area));
    if rf.delay > gz.riseTime
        gz.delay = ceil((rf.delay - gz.riseTime) / opt.system.gradRasterTime) * opt.system.gradRasterTime; % round-up to gradient raster
    end
    if rf.delay < (gz.riseTime + gz.delay)
        rf.delay = gz.riseTime + gz.delay; % these are on the grad raster already which is coarser
    end
end

if rf.ringdownTime > 0
    tFill = (1:round(rf.ringdownTime/1e-6)) * 1e-6;  % Round to microsecond
    rf.t = [rf.t rf.t(end) + tFill];
    rf.signal = [rf.signal, zeros(size(tFill))];
end

rf.signal = rf.signal';


end
