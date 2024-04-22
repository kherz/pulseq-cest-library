import numpy as np
import pypulseq as pp
from write_sequence import write_clinical_sequence_gaussian
import os 


def main():
    # sequence definitions from the bruker sequence
    f_name = 'MRF_CEST_3T_SCONE_007.txt'

    with open(f_name, 'r') as file:
        N = int(file.readline())
        TR, B1, ppm, FA, Tsat = np.zeros((5, N))
        
        for i, line in enumerate(file):
            tr, b1, off, fa, tsat = map(float, line.split())
            TR[i], B1[i], ppm[i], FA[i], Tsat[i] = tr, b1, off, fa, tsat

    # to seconds
    TR = TR / 1000 
    Tsat = Tsat / 1000

    gamma = 267.5153  # [rad / uT]
    seq_fn = 'MRF_CEST_3T_SCONE_007.seq'
    type_s = 'scanner' # type of sequence can be simulator or scanner

    defs = {}

    defs["num_meas"] = N  # number of repetition

    defs["tp"] = 16e-3  # pulse duration [s]
    defs["td"] = 0  # interpulse delay [s]
    defs["dcsat"] = (defs["tp"]) / (defs["tp"] + defs["td"])  # duty cycle

    defs["n_pulses"] = np.ceil(Tsat / (defs["tp"] + defs["td"])).astype(int)  # number of pulses

    defs["tsat"] = Tsat  # saturation time [s]
    defs["trec"] = TR - defs["tsat"]  # net recovery time [s]

    defs["spoiling"] = True

    seqid = os.path.splitext(seq_fn)[1][1:]
    defs['seq_id_string'] = seqid  # unique seq id

    defs['b1pa'] = B1
    defs["b1rms"] = defs["b1pa"]

    defs["offsets_ppm"] = ppm  # offset list [ppm]

    lims = pp.Opts(
        max_grad=40,
        grad_unit="mT/m",
        max_slew=130,
        slew_unit="T/m/s",
        rf_ringdown_time=30e-6,
        rf_dead_time=100e-6,
        rf_raster_time=1e-6,
        gamma=gamma/2/np.pi*1e6,
    )

    defs["gamma_hz"] = lims.gamma * 1e-6
    defs["freq"] = 127.7292
    defs['B0'] = defs['freq'] / defs["gamma_hz"]

    write_clinical_sequence_gaussian(defs, seq_fn, lims)


if __name__ == "__main__":
    # change folder to script directory if needed
    current_directory = os.getcwd()
    script_directory = os.path.dirname(os.path.abspath(__file__))

    if current_directory != script_directory:
        os.chdir(script_directory)

    main()
