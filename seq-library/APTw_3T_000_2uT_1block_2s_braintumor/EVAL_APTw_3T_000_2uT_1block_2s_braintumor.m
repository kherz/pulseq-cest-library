%% EVAL APTw_3T_000_2uT_1block_2s_braintumor
% minimalistic data evaluation for data simulated or acquired with
% APTw_3T_000_2uT_1block_2s_braintumor.seq
%
% Moritz Zaiss 2023

%% read in associated seq file
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;                  % for H [Hz/uT]
seq.read('APTw_3T_000_2uT_1block_2s_braintumor.seq');
defs.offsets_ppm   = seq.definitions('offsets_ppm');
Nmeas=numel(defs.offsets_ppm);

%% read in data from simulation
M_z = load('M_z_APTw_3T_000_2uT_1block_2s_braintumor.seq.txt');

%% read in data from measurement (dicom)
collection = dicomCollection('W:\radiologie\data\MR-Physik\Mitarbeiter\Schuere\Miniskript\dcm\PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015');
V = dicomreadVolume(collection); sz=size(V); V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); V= permute(V,[1 2 4 3]); size(V)
subplot(1,2,1), imagesc(V(:,:,6,1));  subplot(1,2,2), plot(squeeze(V(50,50,1,:)));
sz=size(V);
maskInd=1:sz(1)*sz(2)*sz(3);
V_M_z=permute(V,[4 1 2 3]); M_z=V_M_z(:,maskInd); M_z=double(M_z); figure, plot(M_z(:,5050:5060));

figure, montage1t(squeeze(V_M_z(4,:,:,:)))

%% 
M0=M_z(1,:);
Z=M_z(2:end,:)./M0; 
w=defs.offsets_ppm(2:end);
Z_corr=zeros(size(Z,1),size(Z,2)); dB0_stack=zeros(1,size(Z,2));
% Perform smoothing spline interpolation
tic
for ii=1:size(Z,2)
    try
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
end
toc

figure, plot(w,Z(:,8000),w,Z_corr(:,8000));
% [Z_corrExt_denoised,used_components] = pcca_3D(Z_corrExt,P, (1:numel(P.SEQ.w)), 0, Segment); %principle component analysis (for denoising)
Zref = Z_corr(end:-1:1,:);
MTRasym = Zref-Z_corr;

if size(Z,2)>0
    V_MTRasym=double(V_M_z(2:end,:,:,:))*0;   
    V_MTRasym(:,maskInd)=MTRasym;
    V_Z_corr=double(V_M_z(2:end,:,:,:))*0; 
    V_Z_corr(:,maskInd)=Z_corr; 
end

figure, imagesc(squeeze(V_Z_corr(31,:,:,6)));
figure, imagesc(squeeze(V_MTRasym(31,:,:,6)),[-0.05 0.05]);
w(31)

% colormap(gca,RAINBOW)
