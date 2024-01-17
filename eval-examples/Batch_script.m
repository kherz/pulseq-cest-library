% Batch script examples calling eval scripts

EVAL_APTw_3T('data_flag','re_simulation','bmsim_filename','phantoms\creatine\creatine_3T_pH6.4_T22C_bmsim.yaml','seq_filename', 'APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq');

EVAL_APTw_3T('data_flag','re_simulation','bmsim_filename','phantoms\l-arginin\L-arginin_3T_27mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml','seq_filename', 'APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq');

EVAL_APTw_3T('data_flag','real_data','data_path','W:\radiologie\mr-physik-data\Mitarbeiter\Zaiss_AG\HACKATHON_pulseqCEST_lib\eval\data\PULSEQ_HYBRID_GRE_2_2_5_APTW_001_0009')
