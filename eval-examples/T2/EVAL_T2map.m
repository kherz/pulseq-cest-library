function EVAL_T2map(varargin)

%% EVAL T2map_001_T2prep
% minimalistic data evaluation for T2 prep data acquired with
% T2map_001_T2prep.seq
%
% Moritz Zaiss, Jan-Rüdiger Schüre, Moritz Fabian 
%
% The following flag determines, if you want to operate on real data or
% simulate the data

p = inputParser;
% Define the named parameters and their default values
addParameter(p, 'data_flag', 'real_data');  % simulation, re_simulation or real_data
addParameter(p, 'data_path', '');
addParameter(p, 'bmsim_filename', 'WM_3T_default_7pool_bmsim.yaml');
addParameter(p, 'seq_filename', 'APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq');
parse(p, varargin{:});

data_flag=  p.Results.data_flag;
data_path=  p.Results.data_path;
bmsim_filename=  p.Results.bmsim_filename;
seq_filename=  p.Results.seq_filename;

if strcmp(data_flag,'real_data') && strcmp(data_path,'')
    data_path=uigetdir('','Go to DICOM Directory');
end


%% 1)  read in associated seq file
%get the user specific pulseq path
pulseq_struct=what('pulseq-cest-library');
pulseq_path=pulseq_struct.path;

seq_file_folder_path= [pulseq_path filesep 'seq-library' filesep extractBefore(seq_filename, '.seq')];
seq_file_path= [seq_file_folder_path filesep seq_filename];


% Go to seq file
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;     % for H [Hz/uT]

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

switch data_flag
    case 'simulation'
        %% 2a)  read in data from simulation in pulseq folder
        M_z = load([seq_file_folder_path filesep 'M_z_' seq_filename '.txt']);
    case 're_simulation'
        %% 2b)  re-simulate
        M_z = simulate_pulseqcest(seq_filename, [pulseq_path filesep 'sim-library' filesep 'WM_3T_default_7pool_bmsim.yaml']);
        M_z=M_z';
    case 'real_data'
        %% 2c)  read data from measurement (dicom)
        dcmpath=uigetdir('','Go to DICOM Directory of T2map_001_T2prep  data'); cd(dcmpath)
        collection = dicomCollection(fullfile(dcmpath));
        cd(dcmpath)
        V = double(dicomreadVolume(collection)); sz=size(V); V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); V= permute(V,[1 2 4 3]); size(V)
        
        %Vectorize
        V_M_z=double(permute(V,[4 1 2 3]));       % Changes from 112x92x12x32 to 32x112x92x12
        sz=size(V_M_z);
        maskInd=1:sz(2)*sz(3)*sz(4);
        M_z=V_M_z(:,maskInd);
end

%% 3) get times from seq file
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



%% 4) Evaluation

% Normalization
Z=double(M_z(2:end,:))./double(M_z(1,:));

%reshape segment
Segment=squeeze(V(:, :, :, 1)) > 100;
Segment_resh=reshape(Segment, size(Segment,1)*size(Segment,2)*size(Segment,3),1)';    % size Segment: 112x92x12

% boundaries/start(p0) for:
%    "T2"         "a"      "b"
lb = [  .015      0        -10  ];
ub = [ 1.5          10        10  ];
p0 = [0.05        1         1  ]; %starting values

%define T2 fit function 
t2_fit= @(p,t) p(2)*exp(-t/p(1))+p(3);

%preallocation of Z_fit for speed
Z_fit=zeros(sz(2)*sz(3)*sz(4),1);

for ii = 1:size(Z, 2)
    if all(isfinite(Z(:, ii))) && Segment_resh(ii)==1   
        try
            opts = optimset('Display','off');                        
            [p, ~] = lsqcurvefit(t2_fit,p0,TE',(Z(:,ii)),lb, ub, opts); 
            %disp(p);
            Z_fit(ii,:) = p(1) * 1000; % T2 in ms 
        catch
            disp('something went wrong');
            Z_fit(ii,:) = NaN;
        end
    end
end


%Vectorization Backwards
  sizes=size(V);
  T2map=reshape(Z_fit,[sizes(1) sizes(2) sizes(3)]);

%% 4) Imaging

figure; imagesc(T2map(:,:,6),[0 100]); colormap gray; colorbar; axis image
title('T2map [ms]');

end
