%% EVAL APTw_3T_000_2uT_1block_2s_braintumor
% minimalistic data evaluation for data simulated or acquired with
% APTw_3T_000_2uT_1block_2s_braintumor.seq
%
% Moritz Zaiss 2023

%% 1)  read in associated seq file   
% Go to seq file
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;                  % for H [Hz/uT]
seq.read('APTw_3T_001_2uT_36SincGauss_DC90_2s.seq');
defs.offsets_ppm   = seq.definitions('offsets_ppm');
Nmeas=numel(defs.offsets_ppm);

%% 2a)  read in data from simulation
M_z = load(['M_z_APTw_3T_000_2uT_1block_2s_braintumor.seq.txt']);

%% 2b)  re-simulate
M_z = simulate_pulseqcest('APTw_3T_000_2uT_1block_2s_braintumor.seq','../../sim-library/WM_3T_default_7pool_bmsim.yaml');
M_z=M_z';

%% 2c)  read data from measurement (dicom)
dcmpath=uigetdir('','Go to DICOM Directory'); cd(dcmpath)

% Question if seq file is still the same
question = input('Are the DICOM Files acquired with the same Protocoll Parameters from the PulseqCEST Library? [y/n]','s');
if strcmpi(question, 'y')
    % use above definitions from sction 1)
else
    [seqfile, seqpath]=uigetfile('','*.seq');
    seq.read(fullfile(seqpath,seqfile));
    defs.offsets_ppm   = seq.definitions('offsets_ppm');
    Nmeas=numel(defs.offsets_ppm);
end

cd(dcmpath)
collection = dicomCollection(dcmpath);
V = dicomreadVolume(collection); sz=size(V); V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); V= permute(V,[1 2 4 3]); size(V)
figure;

dcm_image = imagesc(squeeze(V(:, :, 11, 18)));colormap gray
colorbar
hold on;
plot(65,132,'r+')

sizes=size(V);
figure
plot(1:sizes(4),squeeze(V(100,90, 6,:)));

% Vectorization
V_M_z=double(permute(V,[4 1 2 3]));
sz=size(V_M_z);
maskInd=1:sz(2)*sz(3)*sz(4);
M_z=V_M_z(:,maskInd);


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
%Vectorization Backwards
if size(Z,2)>1
    V_MTRasym=double(V_M_z(2:end,:,:,:))*0;
    V_MTRasym(:,maskInd)= MTRasym;
    V_Z_corr=double(V_M_z(2:end,:,:,:))*0;
    V_Z_corr(:,maskInd)= Z_corr;
end


%% 4)  Plots and further graphics
figure;
subplot(1,2,1); plot(w,mean(Z_corr,2),'r.-'); title('Mean Z-spectrum'); set(gca,'Xdir','reverse');
subplot(1,2,2); plot(w,mean(MTRasym,2),'b.-'); title('Mean MTRasym-spectrum'); xlim([0 Inf]);set(gca,'Xdir','reverse');

if size(Z,2)>1
    figure;
    sliceofinterest=6;   % Pick slice for Evaluation
    offsetofinterest=32; % Pick offset for Evaluation
    w(offsetofinterest)
    
    subplot(1,2,1);
    imagesc(squeeze(V_Z_corr(offsetofinterest,:,:,sliceofinterest)),[0.5 1]);  title(sprintf('Z(\\Delta\\omega) = %.2f ppm',w(offsetofinterest)));
    subplot(1,2,2);
    imagesc(squeeze(V_MTRasym(offsetofinterest,:,:,sliceofinterest)),[-0.05 0.05]); title(sprintf('MTRasym(\\Delta\\omega) = %.2f ppm',w(offsetofinterest)));
    %colormap(gca,RAINBOW)
end

