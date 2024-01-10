function EVAL_T1map_001_HypSec

data_flag= 'real_data'; % simulation, re_simulation or real_data
%% 1) Build up Filename, Structure, Paths

%get the user specific pulseq path
pulseq_struct=what('pulseq-cest-library');                  % look in Matlabpath if Folder is already included
pulseq_path=pulseq_struct.path;
seq_filename='T1map_001_3HypSec.seq';
seq_file_folder_path= [pulseq_path filesep 'seq-library' filesep extractBefore(seq_filename, '.seq')];
seq_file_path= [seq_file_folder_path filesep seq_filename];


% Go to seq file
seq = SequenceSBB(getScannerLimits());
gamma_hz  = seq.sys.gamma*1e-6;     % for H [Hz/uT]

question = input('Are the DICOM Files acquired with the same Protocoll Parameters from the PulseqCEST Library? [y/n]','s');
if strcmpi(question, 'y')
   %get the seq file from the library
    seq.read(seq_file_path);
    defs.offsets_ppm = seq.definitions('offsets_ppm');
    Nmeas=numel(defs.offsets_ppm); 
else
    %search your seqfile of the measurement
    [seq_filename, seq_file_folder_path]=uigetfile('','');
    seq.read(fullfile(seq_file_folder_path,seq_filename));
    defs.offsets_ppm = seq.definitions('offsets_ppm');
    Nmeas=numel(defs.offsets_ppm); 
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
        dcmpath=uigetdir('','Go to DICOM Directory'); cd(dcmpath)
   
        collection = dicomCollection(fullfile(dcmpath));
        V= dicomreadVolume(collection); 
        sz=size(V); 

        V=reshape(V,[sz(1) sz(2) Nmeas sz(4)/Nmeas ]); 
        V= permute(V,[1 2 4 3]); size(V)

        % Normalization with TI max
        V_norm=double(V(:,:,:,2:end))./double(V(:,:,:,1));
        size_norm=size(V_norm);

        % Define Segment
        Segment=squeeze(V(:, :, :, 1)) > 100;
        Segment_resh=reshape(Segment, size(Segment,1)*size(Segment,2)*size(Segment,3),1)';    % size Segment: 112x92x12

        % Vectorize
        V_M_z=double(permute(V_norm,[4 1 2 3]));       % Changes from 112x92x12x32 to 32x112x92x12
        sz=size(V_M_z);
        maskInd=1:sz(2)*sz(3)*sz(4);
        M_z=V_M_z(:,maskInd);
end

    %% 3 Evaluation

    % Fitting Parameters 
    P.TI=seq.definitions('TI');P.TI=P.TI.*1000  % conversion from s to ms
    P.B1=seq.definitions('B1rms');
    P.B0=seq.definitions('B0');
    P.Freq=gamma_hz*P.B0;

    P.FIT.options   = [1E-04, 1E-15, 1E-10, 1E-04, 1E-06];
    P.FIT.nIter     = 100;
    P.FIT.modelnum  = 031012;
    P.FIT.extopt=1; % change parameters explicitly
    .FIT.estimatedp0 = 1; %estimate T1 start parameter via linear regression fit

    %     T1          a      c       -->further information in fitmodelfunc_NUM.m
    
    %lb = [0         -80   0     ];
    %ub = [10000      40   40    ];
    p0 = [0.5      1     2000];      % initial Starting values      

    Z_fit = zeros(size(M_z, 2), numel(P.TI)-1);
    t1_stack = zeros(1, size(M_z, 2));

    % Function
    t1_fit = @(p,t)  p(1) - p(2) * exp(-t / p(3));   % Function

    Z=M_z;                                           % Make Copy 
        for ii = 1:size(Z, 2)
            if all(isfinite(Z(:, ii))) && Segment_resh(ii)==1

                try  
                    opts= optimset('Display','off');
                    [p, ~] = lsqcurvefit(t1_fit,p0,P.TI(2:end),(Z(:,ii)),[],[],opts);   
                    %[p, ~] = lsqcurvefit(t1_fit,p0,P.TI(2:end),(Z(:,ii)),lb,ub,opts);  
                    t1_stack(ii) = p(3);    
                catch
                    disp('something went wrong');
                    t1_stack(ii) = NaN;
                end
            end
        end


  %Vectorization Backwards
  sizes=size(V);
  T1=reshape(t1_stack,[sizes(1) sizes(2) sizes(3)])

%% 4 Imaging

 % Display T1
 figure;
 imagesc(T1(:,:,6),[0 2000]);colorbar;title('T1 Map');
end