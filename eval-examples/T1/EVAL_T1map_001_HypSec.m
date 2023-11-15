%% EVAL T1map_001_Hypsec
% minimalistic data evaluation for T1 IR data acquired with
% T1map_001_3HypSec.seq
%
% Moritz Fabian 2023

%% 0) get the seq file infos
pulseq_struct=what('pulseq-cest-library');
pulseq_path=pulseq_struct.path;
seq_filename='T1map_001_3HypSec.seq';
seq_file_folder_path= [pulseq_path filesep 'seq-library' filesep extractBefore(seq_filename, '.seq')];
seq_file_path= [seq_file_folder_path '\' seq_filename];

%initiate seq file
seq = SequenceSBB(getScannerLimits());

question = input('Are the DICOM Files acquired with the same Protocoll Parameters from the PulseqCEST Library? [y/n]','s');
if strcmpi(question, 'y')
   %get the seq file from the library
    seq.read(seq_file_path);
    defs.TI   = seq.definitions('TI');
    Nmeas=numel(defs.TI);
else
    %search your seqfile of the measurement
    [seq_filename, seq_file_folder_path]=uigetfile('','');
    seq.read(fullfile(seq_file_folder_path,seq_filename));
    defs.TI   = seq.definitions('TI');
    Nmeas=numel(defs.TI);
end


%% 1 read data from measurement (dicom)
dcmpath=uigetdir('','Go to DICOM Directory of T1map_001_3HypSec data'); cd(dcmpath)
collection = dicomCollection(fullfile(dcmpath));
cd(dcmpath)
V = double(dicomreadVolume(collection)); sz=size(V); V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); V= permute(V,[1 2 4 3]); size(V)

%% 2) Define Segment
Segment = ones(size(V(:,:,:,1)));
%% 3) get times from seq file
P=logread(seq_filename,seq_file_folder_path);
P.SEQ.w = P.SEQ.TI_list.*1000;% conversion from s to ms as starting values in loadGui are in ms
P.SEQ.FREQ=gamma_*P.SEQ.B0;
P.SEQ.B1=seq.definitions('B1pa');

%% 4) T1 mapping (T1eval_levmar)
TImaxstack=squeeze(V(:,:,:,1));
[image] = NORM_ZSTACK(V,TImaxstack,P,Segment);  %Normalization with max TI
P.SEQ.stack_dim=size(image);

%Fitparameters
P.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
P.FIT.nIter     = 100;
P.FIT.modelnum  = 031012;
P.FIT.extopt=1; % change parameters explicitly
P.FIT.estimatedp0 = 1; %estimate T1 start parameter via linear regression fit
lb = [0         -80   0       ];
ub = [10000      40   40    ];
p0 = [200       0.5     1    ]; %2000 for brain
ROInumber = 1;
P.FIT.lower_limit_fit = lb; P.FIT.upper_limit_fit = ub; P.FIT.start_fit = p0;
 

[popt, P] = FIT_3D(image,P,Segment);   %FIT_3D(Z_stack,P,Segment,slices)
T1map=popt(:,:,:,1);

figure; montage1t(T1map,[500 2000]); colormap gray; colorbar;
title('T1map [ms]');

%% Example to Plot T1map in Nifti File
%nii=load_untouch_nii('MID012001_T1map_001_3HypSec_QUASS.nii')
%nii.img=T1map;
%nii.hdr.dime.dim(5)=1;
%nii.hdr.dime.datatype=16;
%nii.hdr.dime.bitpix=32;
%save_untouch_nii(nii,'T1map');