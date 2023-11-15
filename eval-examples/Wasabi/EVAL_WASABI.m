%% EVAL WASABI 

%% 0 Build up Filename, Structure, Paths

%get the user specific pulseq path
pulseq_struct=what('pulseq-cest-library');                  % look in Matlabpath if Folder is already included
pulseq_path=pulseq_struct.path;
seq_filename='WASABI_3T_001_3p7uT_1block_5ms.seq';
seq_file_folder_path= [pulseq_path filesep 'seq-library' filesep extractBefore(seq_filename, '.seq')];
seq_file_path= [seq_file_folder_path filesep seq_filename];


% Go to seq file
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;     % for H [Hz/uT]
%get the seq file from the library
seq.read(seq_file_path);
defs.offsets_ppm   = seq.definitions('offsets_ppm');
Nmeas=numel(defs.offsets_ppm);  

%% 2a)  read in data from simulation in pulseq folder
M_z = load([seq_file_folder_path filesep 'M_z_' seq_filename '.txt']);

%% 2b)  re-simulate
M_z = simulate_pulseqcest(seq_filename, [pulseq_path filesep 'sim-library' filesep 'WM_3T_default_7pool_bmsim.yaml']);
M_z=M_z';

%% 2c)  read data from measurement (dicom)
dcmpath=uigetdir('','Go to DICOM Directory'); cd(dcmpath)

cd(dcmpath)
collection = dicomCollection(fullfile(dcmpath));

V= dicomreadVolume(collection); 
sz=size(V); 
V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); 
V= permute(V,[1 2 4 3]); size(V)

figure;
dcm_image = imagesc(squeeze(V(:, :, 6, 18)));colormap gray;colorbar

Segment=squeeze(V(:, :, :, 1)) > 100;
Segment_resh=reshape(Segment, size(Segment,1)*size(Segment,2)*size(Segment,3),1)';    % size Segment: 112x92x12

%Vectorize
V_M_z=double(permute(V,[4 1 2 3]));       % Changes from 112x92x12x32 to 32x112x92x12
sz=size(V_M_z);
maskInd=1:sz(2)*sz(3)*sz(4);
M_z=V_M_z(:,maskInd);


%% 3) Evaluation
% Normalization
Z_wasabi=double(M_z(2:end,:))./double(M_z(1,:));


% Fitting Parameters 
B1=seq.definitions('B1');
freq=seq.definitions('FREQ');
w=defs.offsets_ppm(2:end);
t_p=seq.definitions('tp');  
fitoptions=[1E-04, 1E-15, 1E-10, 1E-4, 1E-06];
iterations = 100;
modelnum  = 021031;
c=1;

d=2;


% Compute the analytic function of WASABI curve
wasabi_fit_2abs = @(p,w) abs(p(3) - p(4) .* sin(atan((p(1) ./ (freq / gamma_hz)) ./ (w - p(2)))).^2 .* sin(sqrt((p(1) ./ (freq / gamma_hz)).^2 + (w - p(2)).^2) .* freq .* (2 * pi) .* t_p / 2).^2);
Z=Z_wasabi;

% Adapt the function for use with lsqcurvefit
% This anonymous function captures B1, freq, gamma_, and t_p from the outer scope
%wasabi_fit_function = @(params, offset) wasabi_fit_2abs(params, offset);

% Perform the WASABI fit
dB0_stack = zeros(1, size(Z, 2));
rB1_stack = zeros(1, size(Z, 2));
Z_fit = zeros(size(Z, 2), numel(w));

for ii = 1:size(Z, 2)
    if all(isfinite(Z(:, ii))) && Segment_resh(ii)==1
        try
            p0 = [3.7, -0.1, 1, 2];                              % initial Starting values                 
            %opts = optimset('Display','off');                        
            [p, ~] = lsqcurvefit(wasabi_fit_2abs,p0,w,(Z(:,ii))); 
            disp(p);
            rB1_stack(ii) = p(1) / B1;
            dB0_stack(ii) = p(2);
            Z_fit(ii,:) = wasabi_fit_2abs(p,w);
  
        catch
            disp('something went wrong');
            rB1_stack(ii) = NaN;
            dB0_stack(ii) = NaN;
            Z_fit(ii,:) = NaN;
        end
    end
end


%Vectorization Backwards
  sizes=size(V);
  B1=reshape(rB1_stack,[sizes(1) sizes(2) sizes(3)]);
  B0=reshape(dB0_stack,[sizes(1) sizes(2) sizes(3)]);
  Zfit=reshape(Z_fit,[sizes(1) sizes(2) sizes(3) 31]);
%% 4) plots and further graphics
 % Display B1 and B0
 figure;
 subplot(1,2,1);imagesc(B0(:,:,6),[-0.2 0.2]);colorbar;title('B0 Map'); axis image
 subplot(1,2,2);imagesc(B1(:,:,6),[0.8 1.2]);colorbar;title('B1 Map'); axis image
