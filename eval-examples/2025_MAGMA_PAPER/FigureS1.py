# Figure S1

import sys
from pathlib import Path

original_path = sys.path.copy()
B0_correct = 0

subfigure = input('Please specify which subfigure to generate: ')

if subfigure == 'a' or subfigure == 'b':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("real_data", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.95,
                  (B0_correct-.1))
  sys.path = original_sys_path

elif subfigure == 'c' or subfigure == 'd':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("real_data", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.95,
                  (B0_correct-.05))
  sys.path = original_sys_path

elif subfigure == 'e' or subfigure == 'f':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("real_data", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.95,
                  (B0_correct-.01))
  sys.path = original_sys_path

elif subfigure == 'g' or subfigure == 'h':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("real_data", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.95,
                  B0_correct)
  sys.path = original_sys_path

elif subfigure == 'i' or subfigure == 'j':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("real_data", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.95,
                  (B0_correct+.01))
  sys.path = original_sys_path

elif subfigure == 'k' or subfigure == 'l':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("real_data", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.95,
                  (B0_correct+.05))
  sys.path = original_sys_path

elif subfigure == 'm' or subfigure == 'n':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_APTW_001_RR_0015
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("real_data", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.95,
                  (B0_correct+.1))
  sys.path = original_sys_path
  
else:
  print("Invalid subfigure selected")
