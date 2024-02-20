%% EVAL WASABIT1
% The following flag determines, if you want to operate on real data or
% simulate the data
data_flag= 'real_data'; % simulation, re_simulation or real_data
data_type_flag='nii'; %dcm or nii
%% 1) Build up Filename, Structure, Paths

%get the user specific pulseq path
pulseq_struct=what('pulseq-cest-library');                  % look in Matlabpath if Folder is already included
pulseq_path=pulseq_struct.path;
seq_filename='WASABIT1_3T_001_3p7uT_1block_5p12ms.seq';
seq_file_folder_path= [pulseq_path filesep 'seq-library' filesep extractBefore(seq_filename, '.seq')];
seq_file_path= [seq_file_folder_path filesep seq_filename];


% Go to seq file
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;     % for H [Hz/uT]
%get the seq file from the library
seq.read(seq_file_path);
defs.offsets_ppm   = seq.definitions('offsets_ppm');
Nmeas=numel(defs.offsets_ppm);  

switch data_flag
    case 'simulation'
        %% 2a)  read in data from simulation in pulseq folder
        M_z = load([seq_file_folder_path filesep 'M_z_' seq_filename '.txt']);
    case 're_simulation'
        %% 2b)  re-simulate
        M_z = simulate_pulseqcest(seq_file_path, [pulseq_path filesep 'sim-library' filesep 'WM_3T_default_7pool_bmsim.yaml']);
        M_z=M_z';
       
    case 'real_data'
        switch data_type_flag
            case 'dcm'
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
            case 'nii'
                [niifile, niipath]=uigetfile('*.nii','Get WASABIT1 nifti file.');
                cd(niipath)
                V= readnifti(niifile); 
                sz=size(V); 

                Segment=squeeze(V(:, :, :, 1)) > 100;
                Segment_resh=reshape(Segment, size(Segment,1)*size(Segment,2)*size(Segment,3),1)';    % size Segment: 112x92x12
                
                %Vectorize
                V_M_z=double(permute(V,[4 1 2 3]));       % Changes from 112x92x12x32 to 32x112x92x12
                sz=size(V_M_z);
                maskInd=find(Segment); %new, using the Segment from above
                M_z=V_M_z(:,maskInd);
        end
end

%% 3) Evaluation
% Normalization
Z_wasabi=double(M_z(2:end,:))./double(M_z(1,:));

% Fitting Parameters 
B1pa=seq.definitions('B1pa');
freq=seq.definitions('FREQ');
w=defs.offsets_ppm;
t_p=seq.definitions('tp');
t_rec=seq.definitions('Trec');
fitoptions=[1E-04, 1E-15, 1E-10, 1E-4, 1E-06];
iterations = 100;

%compute w and trec array and name it wt
wt=cat(2,w(2:end)*freq,t_rec(2:end)); %in trec, account for 1st offset being M_0

% Compute the analytic function of WASABIT1 curve
% p(1) = B1, p(2) = dB0, p(3) T1 ;  [3.7, -0.1, 1, 2, 1.5];    
wasabiti_fit_2abs = @(p,w_and_trec_array) abs((1- exp(-w_and_trec_array(:,2)./p(3))) .* (1 - 2 .* (pi * p(1) * gamma_hz * t_p).^2 .* (sinc(t_p .* sqrt((p(1) .* gamma_hz).^2 + (w_and_trec_array(:,1)-p(2)).^2))).^2));


Z=Z_wasabi;

%create a WASABIT1 lookup bib
    B1_bib=[0.8:0.01:1.2]*B1pa;
    B0_bib=[-1:0.05:1]*freq;
    T1_bib=[0.5:0.05:2.5];
    [B1,B0,T1] = ndgrid(B1_bib,B0_bib,T1_bib);
    
    bib_entries=[ B1(:) B0(:) T1(:)];
    wbib=zeros(numel(B1),numel(w(2:end))); %2:end to exclude the M0 offset
    for ii=1:numel(B1)
       wbib(ii,:)= wasabiti_fit_2abs([ B1(ii) B0(ii) T1(ii)], wt);
    end



% Adapt the function for use with lsqcurvefit
% This anonymous function captures B1, freq, gamma_, and t_p from the outer scope
%wasabi_fit_function = @(params, offset) wasabi_fit_2abs(params, offset);

%% Perform the WASABI fit
dB0_stack = zeros(1, size(Z, 2));
rB1_stack = zeros(1, size(Z, 2));
Z_fit = zeros(size(Z, 2), numel(w)-1);
T1_stack = zeros(1, size(Z, 2));

for ii = 1:size(Z, 2)
    if all(isfinite(Z(:, ii)))
        try
            %do a WASABIT1 lookup of start parameters
                SSE=sum((wbib-repmat(Z(:,ii),1,size(wbib,1))').^2,2);
                idx=find(SSE==min(min(min(min(SSE)))),1);
                p0 = bib_entries(idx,:);
                lb=p0-0.5;
                ub=p0+0.5;
                lb(1)=p0(1)*0.2;
                ub(1)=p0(1)*2.5;            
            %p0 = [3.7, -0.1, 1.8];                              % initial Starting values                 
            opts = optimset('Display','off');                        
            [p, ~] = lsqcurvefit(wasabiti_fit_2abs,p0,wt,(Z(:,ii)),lb,ub,opts); 
            rB1_stack(ii) = p(1) / B1pa;
            dB0_stack(ii) = p(2)/freq;
            T1_stack(ii)=p(3);
            Z_fit(ii,:) = wasabiti_fit_2abs(p,wt);
  
        catch
            disp('something went wrong');
            rB1_stack(ii) = NaN;
            dB0_stack(ii) = NaN;
            Z_fit(ii,:) = NaN;
        end
    else
            rB1_stack(ii) = NaN;
            dB0_stack(ii) = NaN;
            Z_fit(ii,:) = NaN; 
    end
end

figure('Name',data_flag);
 plot(mean(Z,2),'rx-', 'DisplayName', 'WASABIT1 measurement'); title('Mean MTRasym-spectrum'); xlim([0 Inf]); hold on;
 plot(mean(Z_fit,1),'k-', 'DisplayName', 'Fit'); title('Mean Z-spectrum'); set(gca,'Xdir','reverse');

if size(Z,2)>1
%Vectorization Backwards
         sizes=size(V);
         B1_map=zeros([sizes(1) sizes(2) sizes(3)]);
         B0_map=zeros([sizes(1) sizes(2) sizes(3)]);
         T1_map=zeros([sizes(1) sizes(2) sizes(3)]);
         Zfit=zeros([sizes(1)*sizes(2)*sizes(3) sizes(4)-1]);
         Z_corrExt=zeros([sizes(1)*sizes(2)*sizes(3) sizes(4)-1]);%for imgui later
      
         B1_map(maskInd)=rB1_stack;
         B0_map(maskInd)=dB0_stack;
         T1_map(maskInd)=T1_stack;
         for ll=1:numel(maskInd)
            Z_corrExt(maskInd(ll),:)=Z_wasabi(:,ll);%for imgui later
            Zfit(maskInd(ll),:)=Z_fit(ll,:);
         end
         Zfit=reshape(Zfit, sizes(1), sizes(2), sizes(3), sizes(4)-1);
         Z_corrExt=reshape(Z_corrExt, sizes(1), sizes(2), sizes(3), sizes(4)-1);

%% 4) Imaging
 % Display B1 and B0
 figure;
 subplot(1,3,1);imagesc(B0_map(:,:,12),[-0.2 0.2]);colorbar;title('B0 Map'); axis image
 subplot(1,3,2);imagesc(B1_map(:,:,12),[0.8 1.2]);colorbar;title('B1 Map'); axis image
 subplot(1,3,3);imagesc(T1_map(:,:,12),[0.8 1.8]);colorbar;title('T1 Map'); axis image


end

%% 5) optional: assign variables to visualize fitting with imgui

%dummy P file
load('WASABI.mat', 'P')

P.SEQ.w=w(2:end);
P.FIT.fitfunc='WASABIT1_FIT_2abs';
P.SEQ.tp=t_p;
P.SEQ.FREQ=freq;
P.SEQ.Trec=t_rec(2:end);
P.SEQ.B1=B1pa;
P.EVAL.w_fit=P.SEQ.w;

Z_uncorr=Z_corrExt; %dummy Z_uncorr
popt(:,:,:,1)=B1_map*P.SEQ.B1;
popt(:,:,:,2)=B0_map;
popt(:,:,:,3)=T1_map;


