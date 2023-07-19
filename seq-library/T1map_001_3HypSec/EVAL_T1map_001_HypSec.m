function EVAL_T1map_001_HypSec

%% 0 Scannerlimits
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;                  % for H [Hz/uT]

%% 1 read data from measurement (dicom)
dcmpath=uigetdir('','Go to DICOM Directory'); cd(dcmpath)
collection = dicomCollection(fullfile(dcmpath));
seriesDescription = "T1map_001_3HypSec_QUASS";                              % Series  Description
indices = find(strcmp(collection.SeriesDescription, seriesDescription));    % Find Dicom File
filePaths = collection.Filenames(indices)
for k=1:numel(filePaths)
V(:,:,:,k) = double(squeeze(dicomreadVolume(filePaths{k,1})));
end

%% Define Segment
Segment = ones(size(V(:,:,:,1)));
%% get times from seq file
seq = SequenceSBB(getScannerLimits());
[seqname seqpath]=uigetfile('','Go to Seq Path','*.seq'); cd(seqpath)
seq.read(seqname)
P=logread(seqname,seqpath)
TI = [10 6 5 4 3 2.5 2 1.5 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1].*1000;  % conversion from s to ms as starting values in loadGui are in ms
P.SEQ.w = TI;
P.SEQ.TI_list = TI;
P.SEQ.FREQ=gamma_hz*P.SEQ.B0;


%% T1 mapping (T1eval_levmar)
TImaxstack=squeeze(V(:,:,:,1))
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