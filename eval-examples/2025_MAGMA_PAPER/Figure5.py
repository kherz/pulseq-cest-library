# Figure 5

import sys
from pathlib import Path

original_path = sys.path.copy()

subfigure = input('Please specify which subfigure to generate: ')

if subfigure == 'a':
  module_path = Path("../../seq-library/WASABI_3T_001_3p7uT_1block_5ms").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from WASABI_3T_001_3p7uT_1block_5ms_parametrized import generate_WASABI_sequence
  generate_WASABI_sequence(11.5,15.5)
  sys.path = original_sys_path

elif subfigure == 'b' or subfigure == 'c' or subfigure == 'd' or subfigure == 'e':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_WASABI_RR_0016
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100

  module_path = Path("../Wasabi").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from WASABI_3T_001_3p7uT_1block_5ms_Eval_parametrized import eval_WASABI
  eval_WASABI("real_data")
  sys.path = original_sys_path
  
elif subfigure == 'f' or subfigure == 'g':
  # Specifications:

  # Data source: PULSEQ_HYBRID_GRE_2_2_5_WASABI_RR_0016
  # 
  # Do you want to specify an ROI [y/n]: y
  # Slice of interest: 5
  # ROI x-minimum: 145
  # ROI x-maximum: 150
  # ROI y-minimum: 95
  # ROI y-maximum: 100
  
  module_path = Path("../Wasabi").resolve()
  original_sys_path = sys.path.copy()
  sys.path.append(str(module_path))
  from WASABI_3T_001_3p7uT_1block_5ms_Eval_parametrized import eval_WASABI
  eval_WASABI("re_simulation")
  sys.path = original_sys_path

else:
  print("Invalid subfigure selected")
