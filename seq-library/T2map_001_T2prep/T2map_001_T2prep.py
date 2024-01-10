# T2map_001_T2prep
# Creates a sequence file for T2 mapping
# T2prep is realized by three pulses of FA 90-180-90 and phases
# 90-180-(-90)
# see Figure 6 of https://doi.org/10.3348/kjr.2017.18.1.113
#
# Tested with pypulseq version 1.3.1post1 and bmctool version 0.6.1
#
# Moritz Zaiss 2023
# moritz.zaiss@fau.de
# --------------------
# Patrick Schuenke 2023
# patrick.schuenke@ptb.de

from pathlib import Path

import numpy as np
import pypulseq as pp
from bmctool.utils.seq.write import write_seq

# get id of generation file
seqid = Path(__file__).stem + "_python"

# get folder of generation file
folder = Path(__file__).parent

# general settings
AUTHOR = "Patrick Schuenke"
FLAG_PLOT_SEQUENCE = False  # plot preparation block?
FLAG_CHECK_TIMING = False  # perform a timing check at the end of the sequence?

# sequence definitions
defs: dict = {}
defs["tp"] = 1.2e-3  # pulse duration [s]
defs["trec"] = 10.0  # recovery time [s]
defs["TE"] = np.array([0, 0.01, 0.025, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0])
defs["offsets_ppm"] = np.zeros(defs["TE"].size)  # offset vector [ppm]
defs["num_meas"] = defs["TE"].size  # number of repetition
defs["b0"] = -1  # dummy B0
defs["b1pa"] = 20  # mean sat pulse B1 [µT]
defs["spoiling"] = "1"  # 0=no spoiling, 1=fully spoiling
defs["seq_id_string"] = seqid  # unique seq id

rb1 = 1.0  # this is just to show the effect of B1 inhomogeneities, e.g. set to rB1=0.8

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

# rf pulses
pre_90 = pp.make_block_pulse(
    flip_angle=rb1 * np.pi / 2,
    duration=defs["tp"],
    system=sys,
)

refocus = pp.make_block_pulse(
    flip_angle=rb1 * np.pi,
    duration=defs["tp"],
    phase_offset=np.pi / 2,
    system=sys,
)

post_90 = pp.make_block_pulse(
    flip_angle=rb1 * np.pi / 2,
    duration=defs["tp"],
    phase_offset=np.pi,
    system=sys,
)

# spoilers
spoil_dur = 5.5e-3  # complete spoiler duration in seconds
spoil_amp = 0.8 * sys.max_grad  # Hz/m
rise_time = 1.0e-3  # spoiler rise time in seconds

gx_spoil, gy_spoil, gz_spoil = [
    pp.make_trapezoid(
        channel=c,
        system=sys,
        amplitude=spoil_amp,
        duration=spoil_dur,
        rise_time=rise_time,
    )
    for c in ["x", "y", "z"]
]

# trec delay
trec_delay = pp.make_delay(defs["trec"])

# pseudo ADC event
pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

# Sequence object
seq = pp.Sequence()

# ===
# RUN
# ===


for m, te in enumerate(defs["TE"]):
    # print progress/offset
    print(f"#{m + 1} / {len(defs['TE'])} : echo time {te:.2f} s")

    real_te_half = te / 2.0 - sys.rf_ringdown_time - sys.rf_dead_time-defs["tp"]

    # add constant trec delay
    if defs["trec"] > 0:
        seq.add_block(trec_delay)

    if te > 0:
        seq.add_block(pre_90)
        seq.add_block(pp.make_delay(real_te_half))
        seq.add_block(refocus)
        seq.add_block(pp.make_delay(real_te_half))
        seq.add_block(post_90)

        seq.add_block(pp.make_delay(100e-6))  # pre-spoil delay
        if defs["spoiling"] == "1" and m > 0:
            seq.add_block(gx_spoil, gy_spoil, gz_spoil)
            seq.add_block(pp.make_delay(100e-6))  # post-spoil delay

    # add pseudo ADC event
    seq.add_block(pseudo_adc)

if FLAG_CHECK_TIMING:
    ok, error_report = seq.check_timing()
    if ok:
        print("\nTiming check passed successfully")
    else:
        print("\nTiming check failed! Error listing follows\n")
        print(error_report)

write_seq(
    seq=seq,
    seq_defs=defs,
    filename=folder / seq_filename,
    author=AUTHOR,
    use_matlab_names=True,
)

# plot the sequence
if FLAG_PLOT_SEQUENCE:
    seq.plot()  # to plot all offsets, remove time_range argument
