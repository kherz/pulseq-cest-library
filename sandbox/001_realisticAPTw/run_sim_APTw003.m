lib_path='D:\root\LABLOG\FAU\MRIlab\SIM\pulseq-cest-library';
seq_path=[lib_path '/seq-library/'];
sim_path=[lib_path '/sim-library/'];

seq_filename=fullfile(seq_path,'/APTw_3T_003_2uT_8block_DC95_834ms_braintumor/APTw_3T_003_2uT_8block_DC95_834ms_braintumor.seq')

seq_filename=fullfile(lib_path,'/sandbox/001_realisticAPTw/APTw_3T_003_2uT_8block_DC95_834ms_braintumor.seq')

% read the .seq-file
seq = mr.Sequence;
seq.read(seq_filename);
% get the definitions in the file
offsets_ppm = seq.definitions('offsets_ppm'); % offsets
m0_offset = seq.definitions('M0_offset');     % m0 offset frequency

  %% call standard sim  for WM and GLIO ( with shifted Lorentzian MT)
figure('Name','WM models');
% M_z = Run_pulseq_cest_Simulation(seq_filename,fullfile(lib_path,'/sim-library/GM_3T_001_bmsim.yaml'));
% Plot_pulseq_cest_Simulation(M_z,offsets_ppm,m0_offset)
% 

file_list= {
    'WM_3T_001_bmsim.yaml'
    'WM_3T_default_7pool_bmsim.yaml'
    'Wang_3T/WM_3T_Wang2020_5pool_bmsim.yaml'
    'Stanisz_3T//WM_3T_Stanisz2005_5pool_bmsim.yaml'  
    'vanZijl2018_3T/WM_3T_003_bmsim.yaml'
    'Heo2019_3T/WM_3T_Heo2019_4pool_bmsim.yaml'   
    'Liu2013_7T/WM_3T_Liu2013_4pool_bmsim.yaml'
    'Glang2022/BF_volTueb1_3T_Lorentz_wo_0ppm_4cestpools_WM.yaml'      
    };

shortID= {
    'WM_3T_001'
    'WM_3T_default_7pool'
    'Wang_3T'
    'Stanisz_3T'
    'vanZijl2018_3T'
    'Heo2019_3T'    
    'Liu2013_7T'
    'Glang2022'    
    };

for ii=1:length(file_list)
    bmsimfilename=file_list{ii};
M_z = simulate_pulseqcest(seq_filename,[sim_path bmsimfilename]);
plotSimulationResults(M_z,offsets_ppm,m0_offset);
end

legend(shortID, 'Interpreter', 'none');

%%

uiopen('D:\root\LABLOG\FAU\11_CESTpulseq_standard\APTw_standard\multi-site-measurment_2\Z-spec_smoothspline_and_linear.fig',1)

