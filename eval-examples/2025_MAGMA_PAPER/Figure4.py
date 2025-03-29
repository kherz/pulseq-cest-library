# Figure 4 

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

elif subfigure == 'b' or subfigure == 'c' or subfigure == 'd':
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

elif subfigure == 'e' or subfigure == 'f' or subfigure == 'g':
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("re_simulation", 
                  "",
                  "phantoms/creatine/creatine_3T_pH6.4_T22C_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.99,
                  0)
  sys.path = original_sys_path

elif subfigure == 'h' or subfigure == 'i' or subfigure == 'j':
  module_path = Path("../MTRasym").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from APTw_3T_eval_parametrized import evaluate_APTw_3T
  evaluate_APTw_3T("re_simulation", 
                  "",
                  "WM_3T_001_bmsim.yaml",
                  "APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq",
                  "spline",
                  0.99,
                  0)
  sys.path = original_sys_path

else:
  print("Invalid subfigure selected")
