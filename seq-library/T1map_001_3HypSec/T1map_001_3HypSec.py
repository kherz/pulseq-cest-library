# T1map_001_3HypSec
# T1 Saturation Recovery Protocol with 3 adiabatic HypSec pulses
#
# IMPORTANT INFO: if you get an AssertionError in the timing check, just comment out the corresponding assert statement:
#       "assert abs(duration - block_duration) < eps"
# in pypulseq\Sequence\block.py. This is a known bug in pypulseq 1.3.1post1.
#
# Tested with pypulseq version 1.3.1post1 and bmctool version 0.6.0
#
# Patrick Schuenke 2023
# patrick.schuenke@ptb.de

from pathlib import Path

import numpy as np
import pypulseq as pp
from bmctool.utils.pulses.make_hypsec_half_passage import make_hypsec_half_passage_rf
from bmctool.utils.seq.write import write_seq

# get id of generation file
seqid = Path(__file__).stem + "_python"

# get folder of generation file
folder = Path(__file__).parent

# general settings
AUTHOR = "Patrick Schuenke"
FLAG_PLOT_SEQUENCE = False  # plot preparation block?
FLAG_CHECK_TIMING = False  # perform a timing check at the end of the sequence?
FLAG_POST_PREP_SPOIL = True  # add spoiler after preparation block?

# sequence definitions
defs: dict = {}
defs["b1pa"] = 20  # B1 peak amplitude [µT] (b1rms calculated below)
defs["b1rms"] = defs["b1pa"]
defs["b0"] = 1  # B0 [T]
defs["n_pulses"] = 3  # number of pulses  #
defs["tp"] = 8e-3  # pulse duration [s]
defs["trec"] = 1.0  # recovery time [s]
defs["TI"] = np.array([10, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])

defs["offsets_ppm"] = np.zeros(defs["TI"].size)  # offset vector [ppm]
defs["num_meas"] = defs["offsets_ppm"].size  # number of repetition
defs["seq_id_string"] = seqid  # unique seq id
defs["spoiling"] = "1" if FLAG_POST_PREP_SPOIL else "0"

seq_filename = defs["seq_id_string"] + ".seq"

# scanner limits
sys = pp.Opts(
    max_grad=40,
    grad_unit="mT/m",
    max_slew=130,
    slew_unit="T/m/s",
    rf_ringdown_time=30e-6,
    rf_dead_time=100e-6,
    rf_raster_time=1e-6,
    gamma=42576400,
)

# ===========
# PREPARATION
# ===========

# spoilers
spoil_amp0 = 0.8 * sys.max_grad  # Hz/m
spoil_amp1 = -0.7 * sys.max_grad  # Hz/m
spoil_amp2 = 0.6 * sys.max_grad  # Hz/m

rise_time = 1.0e-3  # spoiler rise time in seconds
spoil_dur = 5.5e-3  # complete spoiler duration in seconds

gx_spoil0, gy_spoil0, gz_spoil0 = [
    pp.make_trapezoid(channel=c, system=sys, amplitude=spoil_amp0, duration=spoil_dur, rise_time=rise_time)
    for c in ["x", "y", "z"]
]
gx_spoil1, gy_spoil1, gz_spoil1 = [
    pp.make_trapezoid(channel=c, system=sys, amplitude=spoil_amp1, duration=spoil_dur, rise_time=rise_time)
    for c in ["x", "y", "z"]
]
gx_spoil2, gy_spoil2, gz_spoil2 = [
    pp.make_trapezoid(channel=c, system=sys, amplitude=spoil_amp2, duration=spoil_dur, rise_time=rise_time)
    for c in ["x", "y", "z"]
]

# RF pulses
hs_pulse = make_hypsec_half_passage_rf(amp=defs["b1rms"], system=sys)

# pseudo ADC event
pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

# delays
trec_delay = pp.make_delay(defs["trec"])

# Sequence object
seq = pp.Sequence()

# ===
# RUN
# ===


for m, t_inv in enumerate(defs["TI"]):
    # print progress/offset
    print(f"#{m + 1} / {len(defs['TI'])} : inversion time {t_inv:.2f} s")

    # reset accumulated phase
    accum_phase = 0

    # add constant trec delay
    if defs["trec"] > 0:
        seq.add_block(trec_delay)

    # add preparation (adiabatic excitation + spoiler) block
    for i in range(defs["n_pulses"]):
        seq.add_block(hs_pulse)
        if FLAG_POST_PREP_SPOIL:
            if i % 3 == 0:
                seq.add_block(gx_spoil0, gy_spoil1, gz_spoil2)
            elif i % 2 == 0:
                seq.add_block(gx_spoil2, gy_spoil0, gz_spoil1)
            else:
                seq.add_block(gx_spoil1, gy_spoil2, gz_spoil0)

    # add variable inversion time delay
    seq.add_block(pp.make_delay(t_inv))

    # add pseudo ADC event
    seq.add_block(pseudo_adc)

if FLAG_CHECK_TIMING:
    ok, error_report = seq.check_timing()
    if ok:
        print("\nTiming check passed successfully")
    else:
        print("\nTiming check failed! Error listing follows\n")
        print(error_report)

write_seq(seq=seq, seq_defs=defs, filename=folder / seq_filename, author=AUTHOR, use_matlab_names=True)

# plot the sequence
if FLAG_PLOT_SEQUENCE:
    seq.plot()  # to plot all offsets, remove time_range argument
