# Figure 2 

import sys
from pathlib import Path

original_path = sys.path.copy()

subfigure = input('Please specify which subfigure to generate: ')

if subfigure == 'a':
  module_path = Path("../../seq-library/APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor_parametrized import generate_APT_sequence
  generate_APT_sequence(3,6)
  sys.path = original_sys_path

elif subfigure == 'b' or subfigure == 'c' or subfigure == 'd' or subfigure == 'f' or subfigure == 'g':
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
                  0)
  sys.path = original_sys_path

elif subfigure == 'e':
  print("In order to generate this specialized figure, please edit the following lines in APTw_3t_eval_parametrized:\n \
        Z = m_z / M0 --> Z = m_z \n \
        plt.imshow(V_Z_corr_reshaped[:, :, slice_of_interest, offset_of_interest], vmin=0.5, vmax=1) --> plt.imshow(V_Z_corr_reshaped[:, :, slice_of_interest, offset_of_interest])")
  print("Then run this code again and specify subfigure f (since e is a specialized version of subfigure f)")

elif subfigure == 'h' or subfigure == 'i' or subfigure == 'j':
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
  evaluate_APTw_3T("re_simulation", 
                  "",
                  "phantoms/l-arginin/L-arginin_3T_20mM_pH4_T1_1500ms_T2_1000ms_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.99,
                  0)
  sys.path = original_sys_path
  
else:
  print("Invalid subfigure selected")
