function EVAL_APTw_3T(varargin)
%% EVAL APTw_3T_000_2uT_1block_2s_braintumor
% minimalistic data evaluation for data simulated or acquired with
% APTw_3T_000_2uT_1block_2s_braintumor.seq
% Inputs:
%   data_flag - A string parameter to specify the type of data to process.
%   it determines if you want to operate on real data or simulate the data
%       Options are:
%       'simulation'     - for simulated data as provided in the pulseq-cest seq folder
%       're_simulation'  - for re-simulated data using the bmsim.yaml file
%       'real_data'      - for real data. (Default)
%
%   data_path - A string specifying the path to the data.
%               If left empty, a dialog box will prompt for a directory.
%
p= inputParser;
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
%get the seq file from the library
disp('Read seq-file ...');
seq.read(seq_file_path);
defs.offsets_ppm   = seq.definitions('offsets_ppm');
Nmeas=numel(defs.offsets_ppm);

switch data_flag
    case 'simulation'
        %% 2a)  read in data from simulation in pulseq folder
        disp('Evaluating existing simulation data of provided Pulseq-CEST library seq-file');
        M_z = load([seq_file_folder_path filesep 'M_z_' seq_filename '.txt']);
        
    case 're_simulation'
        %% 2b)  re-simulate
        disp('Simulating and Evaluating simulation data using seq and yaml file.');
        M_z = simulate_pulseqcest([seq_file_folder_path filesep seq_filename], [pulseq_path filesep 'sim-library' filesep bmsim_filename]);
        M_z=M_z';
        
    case 'real_data'
        %% 2c)  read data from measurement (dicom)
        % Question if seq file is still the same
        disp('Evaluating real data... It is assumed that the DICOM Files were acquired with the standrad or provided seq file.');
        collection = dicomCollection(data_path);
        V = dicomreadVolume(collection); sz=size(V); V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); V= permute(V,[1 2 4 3]); size(V)
        
        % Vectorization
        V_M_z=double(permute(V,[4 1 2 3]));
        sz=size(V_M_z);
        maskInd=1:sz(2)*sz(3)*sz(4);
        M_z=V_M_z(:,maskInd);
end

%% 3) Evaluation
M0=M_z(1,:);
M0(M0<0.2*mean(M0,2))=0;  % filter pixels that gave less than 10% of the mean intensity
Z=M_z(2:end,:)./M0; % Normalization
w=defs.offsets_ppm(2:end);
Z_corr=zeros(size(Z,1),size(Z,2)); dB0_stack=zeros(1,size(Z,2));
% % Perform smoothing spline interpolation
for ii=1:size(Z,2)    % min Z, B0 Correction
    if  isfinite(Z(:,ii))
        pp = csaps(w, Z(:,ii), 0.95);
        w_fine=-1:0.005:1;
        Z_fine = ppval(pp, w_fine);
        [~, MINidx] = min(Z_fine);
        dB0=w_fine(MINidx);
        dB0_stack(1,ii)=dB0;
        Z_corr(:,ii)= ppval(pp, w+dB0);
        %         disp(ii/size(Z,2));
    end
end

% Denoising here
% [st, rs] = system('git clone https://github.com/cest-sources/CEST-AdaptiveDenoising --depth 1');
% AdaptiveDenoising(Data,Segment,CriterionOrNumber,verbosity);
% [Z_corrExt_denoised,used_components] = pcca_3D(Z_corrExt,P, (1:numel(P.SEQ.w)), 0, Segment); %principle component analysis (for denoising)

%  Calc Zspec and MTRasym
Zref = Z_corr(end:-1:1,:);
MTRasym = Zref-Z_corr;
FS_MTRasym = MTRasym.*0.5^2./(Zref.^2);
%Vectorization Backwards
if size(Z,2)>1
    V_MTRasym=double(V_M_z(2:end,:,:,:))*0;
    V_MTRasym(:,maskInd)= MTRasym;
    V_Z_corr=double(V_M_z(2:end,:,:,:))*0;
    V_Z_corr(:,maskInd)= Z_corr;
    V_FS_MTRasym=double(V_M_z(2:end,:,:,:))*0;
    V_FS_MTRasym(:,maskInd)= FS_MTRasym;
end


%% 4)  Visualization
figure('Name',data_flag);
subplot(1,2,1); plot(w,mean(Z_corr,2),'r.-', 'DisplayName', 'Measurement'); title('Mean Z-spectrum'); set(gca,'Xdir','reverse');
subplot(1,2,2); plot(w,mean(MTRasym,2),'r.-', 'DisplayName', 'Measurement'); title('Mean MTRasym-spectrum'); xlim([0 Inf]);set(gca,'Xdir','reverse');

if size(Z,2)>1
    figure;
    sliceofinterest=6;   % Pick slice for Evaluation
    offsetofinterest=31; % Pick offset for Evaluation
    w(offsetofinterest)
    
    subplot(2,2,1);
    imagesc(squeeze(V_Z_corr(offsetofinterest,:,:,sliceofinterest)),[0.5 1]);  title(sprintf('Z(\\Delta\\omega) = %.2f ppm',w(offsetofinterest)));colorbar;
    subplot(2,2,2);
    imagesc(squeeze(V_MTRasym(offsetofinterest,:,:,sliceofinterest)),[-0.05 0.05]); title(sprintf('MTRasym(\\Delta\\omega) = %.2f ppm',w(offsetofinterest)));colorbar;
    subplot(2,2,4);
    imagesc(squeeze(V_FS_MTRasym(offsetofinterest,:,:,sliceofinterest)),[-0.05 0.05]); title(sprintf('FS-MTRasym(\\Delta\\omega) = %.2f ppm',w(offsetofinterest)));colorbar;
    
    %colormap(gca,RAINBOW)
end

end

