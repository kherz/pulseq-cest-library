%% T2 prep sequence

% Adding path to 3rdParty folder with pulseq master
addpath(genpath('../3rdParty'))

% specify which sequence you want to create
sequenceType = 'SCANNER_PULSEQ'; % SCANNER or SIMULATION

tPrepTimes = [0 0.01 0.025 0.03 0.04 0.05 0.1 0.2 0.3 0.5 1.0];
%tPrepTimes = [0.01 0.025 0.03];
numPreps = numel(tPrepTimes);
tRec = 10;

dummyPulse = false;

outputFileName = ['T2prep_' num2str(numPreps) '_9_11.seq'];
%%
if strcmp(sequenceType, 'SCANNER_PULSEQ')
    externalSeqFileName = 'seqs/SCANNER_PULSEQ.seq';
    if ~exist(externalSeqFileName, 'file')
        error([externalSeqFileName ' not found!']);
    end
else
    error(['this works only for the pulseq sbb']);
end

%% pulse parameters

gyroRatio = 42.5764*2*pi;
fieldStrength = 3; % T

% set the spoiling for the scanner sequence
spoiling = 1; % 0=no spoiling, 1=before readout

% Set system limits (Siemens Prisma: MaxGrad=80mT/m, MaxSlew=200T/m/s (?))
lims = mr.opts('MaxGrad',40,'GradUnit','mT/m',...
    'MaxSlew',130,'SlewUnit','T/m/s', ...
    'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6, 'rfRasterTime',1e-6);
% SL pulses are played out during the interpulse delay. We need to make

%% load readout sequence
externalSequence = mr.Sequence();
externalSequence.read(externalSeqFileName);
readoutTime = externalSequence.duration;

%% spoiler objects
seq=mr.Sequence();              % Create a new sequence object

% spoilers
spoilAmplitude = 0.8 .* lims.maxGrad; % in Hz/m % FG: hard coded to 80% of maximum gradient, more should be possible
riseTime = 1e-3;
spoilDuration = 4500e-6+riseTime; % in s
gxSpoil=mr.makeTrapezoid('x','Amplitude',spoilAmplitude,'Duration',spoilDuration,'riseTime', riseTime,'system',lims);
gySpoil=mr.makeTrapezoid('y','Amplitude',spoilAmplitude,'Duration',spoilDuration,'riseTime', riseTime,'system',lims);
gzSpoil=mr.makeTrapezoid('z','Amplitude',spoilAmplitude,'Duration',spoilDuration,'riseTime', riseTime,'system',lims);

%% build the sequence
PulseDuration =0.0012; % [s] % Single saturation pulse length
preEchoPulse = mr.makeBlockPulse(pi/2, 'Duration', PulseDuration, 'system',lims);
echoPulse = mr.makeBlockPulse(pi, 'Duration', PulseDuration,'Phase', pi/2,'system',lims);
postEchoPulse = mr.makeBlockPulse(pi/2, 'Duration', PulseDuration,'Phase', pi,'system',lims);

if dummyPulse
   % seq.addBlock(mr.makeBlockPulse(deg2rad(90),'Duration', 1e-3, 'freqOffset', -4000, 'system',lims));
    seq.addBlock(preEchoPulse);
end

for ii =1:numPreps
    realTehalf = tPrepTimes(ii)/2-lims.rfRingdownTime-lims.rfDeadTime;
    % recover time
    
    seq.addBlock(mr.makeDelay(tRec));
    if tPrepTimes(ii) > 0 % m0 tthe beginning
        
        seq.addBlock(preEchoPulse);
        seq.addBlock(mr.makeDelay(realTehalf));
        seq.addBlock(echoPulse);
        seq.addBlock(mr.makeDelay(realTehalf));
        seq.addBlock(postEchoPulse);
        
        seq.addBlock(mr.makeDelay(100e-6)); % prespoiler
        if spoiling == 1 && ii > 1 % spoil before readout after saturation
            seq.addBlock(gxSpoil,gySpoil,gzSpoil);
            seq.addBlock(mr.makeDelay(100e-6)); % necessary?
        end
    end
    seq.addBlock(mr.makeAdc(1, 'Duration', 1e-3));
%     %%% Imaging block, take from input seq file
%     for kk=1:length(externalSequence.blockEvents)
%         seq.addBlock(externalSequence.getBlock(kk));
%     end
    %%%%%%%%%%%%%%%%%%%%%
    
end

%%
seq.write(outputFileName);   % Output sequence for scanner

if ~seq.isvalid
    warning ('Sequence did not pass validation check!');
end

%% Ploting the resulting seq
% Temp_Seq_for_ploting = mr.Sequence();
% Temp_Seq_for_ploting.read('MRF_CEST_SCANNER_51_Amide_with_M0_3T.seq');
% Temp_Seq_for_ploting.plot

