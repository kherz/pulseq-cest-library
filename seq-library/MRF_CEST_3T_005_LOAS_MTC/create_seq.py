from sequences import create_defs_loas, write_sequence_clinical
import os
import numpy as np


def main():
    lens = [10, 40]
    seq_fn = 'MRF_CEST_3T_LOAS_MTC_005'
    fa = 90 

    for l in lens:
        defs, lims = create_defs_loas(l, clinical=True)
        fn = f'{seq_fn}_{l}.seq'
        TR = defs['tsat'] + defs['trec']
        defs['seq_id_string'] = f'{seq_fn}_{l}'
        write_sequence_clinical(defs, fn, lims)


if __name__ == '__main__':
    # change folder to script directory if needed
    current_directory = os.getcwd()
    script_directory = os.path.dirname(os.path.abspath(__file__))

    if current_directory != script_directory:
        os.chdir(script_directory)

    main()
