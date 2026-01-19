# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 17:58:21 2025

@author: schuerjn
"""

## Go to eval:_example Folder. Make sure you are in there.

from pathlib import Path
import yaml
import os
from io import StringIO

from bmctool.parameters import Parameters
from MTRasym.EVAL_APTw_3T import EVAL_APTw_3T
from Wasabi.EVAL_WASABI import EVAL_WASABI
from T1.Eval_T1map_001_HypSec import EVAL_T1
from T2.Eval_T2map_001_HypSec import EVAL_T2





# Funktion: YAML-Datei laden
def load_yaml(file_path):
    with open(file_path, 'r') as file:
        return yaml.safe_load(file)
   

if __name__ == "__main__":
    # Original-YAML-Datei laden
    evalname='MTRasym'
    bmsim_filename = 'WM_3T_default_7pool_bmsim.yaml'
    yaml_name = Path(bmsim_filename)   
    yaml_path = Path.cwd().parent/ "sim-library" / bmsim_filename
    
    # load config file(s)
    yaml_struc = Parameters.from_yaml(yaml_path)
    
    
    # Change Parameters in yaml file  
    yaml_struc.system.b0_inhom = 0  # Beispieländerung
   
    # Changes in P struc
    #P = {
    #"smooth_param": 0.95
    #}
    
    
    script_path = Path.cwd().parent/ "eval-examples" / evalname
    os.chdir(script_path)
    # Skript aufrufen mit dem modifizierten YAML-Inhalt
    EVAL_APTw_3T(
        data_flag='real_data',
        data_path='', 
        bmsim_filename= yaml_struc,  # Übergebe den StringIO-Objekt als Dateiinhalt   
        seq_filename='APTw_3T_001_2uT_36SincGauss_DC90_2s_braintumor.seq',
        P_Eval_smoothness=0.95,
    )
    
   
#    EVAL_WASABI(
#        data_flag='re_simulation',
#        data_path='', 
#        bmsim_filename= 'WM_3T_default_7pool_bmsim.yaml',  # Übergebe den StringIO-Objekt als Dateiinhalt   
#        seq_filename='WASABI_3T_001_3p7uT_1block_5ms.seq',
#        P_Eval_smoothness=0.95
#    )
#    
#    
#    EVAL_T1(
#        data_flag='re_simulation',
#        data_path='', 
#        bmsim_filename= 'WM_3T_default_7pool_bmsim.yaml',  # Übergebe den StringIO-Objekt als Dateiinhalt   
#        seq_filename='T1map_001_3HypSec.seq'
#    )
#    
#    
#    EVAL_T2(
#        data_flag='re_simulation',
#        data_path='', 
#        bmsim_filename= 'WM_3T_default_7pool_bmsim.yaml',  # Übergebe den StringIO-Objekt als Dateiinhalt   
#        seq_filename='T2map_001_T2prep.seq'
#    )