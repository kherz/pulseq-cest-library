%% EVAL T2map_001_T2prep
% minimalistic data evaluation for T2 prep data acquired with
% T2map_001_T2prep.seq
%
% Moritz Fabian 2023

%% 0) get the seq file infos
pulseq_struct=what('pulseq-cest-library');
pulseq_path=pulseq_struct.path;
seq_filename='T2map_001_T2prep.seq';
seq_file_folder_path= [pulseq_path '\seq-library\' extractBefore(seq_filename, '.seq')];
seq_file_path= [seq_file_folder_path '\' seq_filename];

%initiate seq file
seq = SequenceSBB(getScannerLimits());

question = input('Are the DICOM Files acquired with the same Protocoll Parameters from the PulseqCEST Library? [y/n]','s');
if strcmpi(question, 'y')
   %get the seq file from the library
    seq.read(seq_file_path);
    defs.TE   = seq.definitions('TE');
    Nmeas=numel(defs.TE);
else
    %search your seqfile of the measurement
    [seq_filename, seq_file_folder_path]=uigetfile('','');
    seq.read(fullfile(seq_file_folder_path,seq_filename));
    defs.TE   = seq.definitions('TE');
    Nmeas=numel(defs.TE);
end

%% 1 read data from measurement (dicom)
dcmpath=uigetdir('','Go to DICOM Directory of T2map_001_T2prep  data'); cd(dcmpath)
collection = dicomCollection(fullfile(dcmpath));
cd(dcmpath)
V = double(dicomreadVolume(collection)); sz=size(V); V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); V= permute(V,[1 2 4 3]); size(V)
%% 2 Define Segment
Segment = ones(size(V(:,:,:,1)));
%% 3) get times from seq file
P_T2=logread(seq_filename,seq_file_folder_path);

TE = [];
num_adc = 0;
% loop through seq and add times to get exact TE value
nB = 1;
while nB < numel(seq.blockEvents)
    block = seq.getBlock(nB);
    if ~isempty(block.rf)
        num_adc = num_adc + 1;
        % 1st block: from middle of 1st RF pulse to end of pulse
        % including delay for hardware at the end
        c_te = numel(find(block.rf.signal>1e-6))/2*1e-6; % middle of pulse
        c_te = c_te + numel(find(block.rf.signal<1e-6))*1e-6; % hardware related delay
        nB = nB+1;
        % 2nd block: delay between 90 and 180 pulse
        block = seq.getBlock(nB);
        c_te = c_te + block.delay.delay;
        nB = nB+1;
        % 3rd block: 180 pulse including its hardware related delay
        block = seq.getBlock(nB);
        c_te = c_te + block.rf.delay + numel(block.rf.signal)*1e-6;
        nB = nB+1;
        % 4th block: delay between 180 and -90 pulse
        block = seq.getBlock(nB);
        c_te = c_te + block.delay.delay;
        nB = nB+1;
        % 5th block: hardware related delay of -90 pulse + half the duration
        % of the pulse
        block = seq.getBlock(nB);
        c_te = c_te + block.rf.delay + numel(find(block.rf.signal>1e-6))/2*1e-6;
        TE(num_adc)= c_te;
    end
    nB = nB+1;
end

%remove first readout in T2 stack
image = double(V);
T2_stack = image(:,:,:,2:end)./image(:,:,:,1);
P_T2.SEQ.w = TE;

%% just for now .. (needed in levmar_fit.m) / may be removed later on
P_T2.SEQ.FREQ = [];
P_T2.SEQ.tp = [];
P_T2.SEQ.B1 = [];

% information about fit
P_T2.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
P_T2.FIT.nIter     = 200;
P_T2.FIT.modelnum  = 051011; % multi echo T2 function: 
P_T2.FIT.extopt=1;   % change parameters explicitly

% boundaries/start(p0) for:
%    "T2"         "a"      "b"
lb = [  .015      0        -10  ];
ub = [ 1          10        10  ];
p0 = [0.05        1         1  ];
P_T2.FIT.lower_limit_fit = lb; P_T2.FIT.upper_limit_fit = ub; P_T2.FIT.start_fit = p0;

tic ;
popt= FIT_3D(T2_stack,P_T2,Segment);
T2map=popt(:,:,:,1);
toc
T2map=T2map.*1000;  % Converting from from s in ms

figure; montage1t(T2map,[0 1000]); colormap gray; colorbar;
title('T2map [ms]');



