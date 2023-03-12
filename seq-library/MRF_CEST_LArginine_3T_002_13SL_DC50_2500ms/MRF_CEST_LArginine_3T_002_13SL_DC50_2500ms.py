# MRF_CEST_LArginine_3T_002_13SL_DC50_2500ms
# L-Arginine CEST Fingerprinting protocol
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
defs["b1pa"] = np.array(
    [
        2.0,
        2.0,
        1.7,
        1.5,
        1.2,
        1.2,
        3.0,
        0.5,
        3.0,
        1.0,
        2.2,
        3.2,
        1.5,
        0.7,
        1.5,
        2.2,
        2.5,
        1.2,
        3.0,
        0.2,
        1.5,
        2.5,
        0.7,
        4.0,
        3.2,
        3.5,
        1.5,
        2.7,
        0.7,
        0.5,
    ]
)
defs["b1rms"] = defs["b1pa"]
defs["b0"] = 3  # B0 [T]
defs["n_pulses"] = 13  # number of pulses
defs["num_meas"] = 30  # number of repetition
defs["tp"] = 100e-3  # pulse duration [s]
defs["td"] = 100e-3  # interpulse delay [s]
defs["trec"] = 1.0  # recovery time [s]
defs["offsets_ppm"] = np.ones(defs["num_meas"]) * 3.0

defs["dcsat"] = (defs["tp"]) / (defs["tp"] + defs["td"])  # duty cycle

defs["tsat"] = defs["n_pulses"] * (defs["tp"] + defs["td"]) - defs["td"]  # saturation time [s]
defs["seq_id_string"] = seqid  # unique seq id
defs["spoiling"] = "1" if FLAG_POST_PREP_SPOIL else "0"

seq_filename = defs["seq_id_string"] + ".seq"

# scanner limits
sys = pp.Opts(
    max_grad=40,
    grad_unit="mT/m",
    max_slew=130,
    slew_unit="T/m/s",
    rf_ringdown_time=50e-6,
    rf_dead_time=200e-6,
    rf_raster_time=1e-6,
    gamma=42576400,
)

GAMMA_HZ = sys.gamma * 1e-6
defs["freq"] = defs["b0"] * GAMMA_HZ  # Larmor frequency [Hz]

# ===========
# PREPARATION
# ===========

# spin lock specific preparation
# td should be between block pulses, so we fit the tipping pulses in the
# delay between pulses
tp_sl = 1e-3  # duration of tipping pulse for sl
td_sl = sys.rf_dead_time + sys.rf_ringdown_time  # delay between tip and sat pulse
sl_time_per_sat = 2 * (tp_sl + td_sl)  # additional time of sl pulses for 1 sat pulse
assert (defs["trec"] >= sl_time_per_sat, "DC too high for SL prepatration pulses!")


# spoiler
spoil_amp = 0.8 * sys.max_grad  # Hz/m
rise_time = 1.0e-3  # spoiler rise time in seconds
spoil_dur = 6.5e-3  # complete spoiler duration in seconds

gx_spoil, gy_spoil, gz_spoil = [
    pp.make_trapezoid(channel=c, system=sys, amplitude=spoil_amp, duration=spoil_dur, rise_time=rise_time)
    for c in ["x", "y", "z"]
]

# pseudo ADC event
pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

# delays
trec_delay = pp.make_delay(defs["trec"])

# Sequence object
seq = pp.Sequence()

# ===
# RUN
# ===

offsets_hz = defs["offsets_ppm"] * defs["freq"]  # convert from ppm to Hz

for m, b1 in enumerate(defs["b1rms"]):
    # print progress/offset
    print(f"#{m + 1} / {len(offsets_hz)} : offset {offsets_hz[m] / defs['freq']:.2f} ppm ({offsets_hz[m]:.3f} Hz)")

    # reset accumulated phase
    accum_phase = 0

    # add delay
    if defs["trec"] > 0:
        seq.add_block(trec_delay)

    # prep and set rf pulse
    flip_angle_sat = defs["b1pa"][m] * GAMMA_HZ * 2 * np.pi * defs["tp"]
    sat_pulse = pp.make_block_pulse(flip_angle=flip_angle_sat, duration=defs["tp"], system=sys)
    sat_pulse.freq_offset = offsets_hz[m]
    accum_phase = np.mod(offsets_hz[m] * 2 * np.pi * defs["tp"], 2 * np.pi)

    # prep tipping pulses
    flip_angle_tip = np.arctan(b1 / (defs["offsets_ppm"][m] * defs["b0"]))
    pre_sl_pulse = pp.make_block_pulse(flip_angle=flip_angle_tip, duration=tp_sl, phase_offset=-(np.pi / 2), system=sys)
    post_sl_pulse = pp.make_block_pulse(
        flip_angle=flip_angle_tip, duration=tp_sl, phase_offset=accum_phase + (np.pi / 2), system=sys
    )
    # sl phase cycling
    phase_cycling = 50 / 180 * np.pi

    for n in range(defs["n_pulses"]):
        pre_sl_pulse.phase_offset = pre_sl_pulse.phase_offset + phase_cycling
        sat_pulse.phase_offset = sat_pulse.phase_offset + phase_cycling
        post_sl_pulse.phase_offset = post_sl_pulse.phase_offset + phase_cycling
        seq.add_block(pre_sl_pulse)
        seq.add_block(sat_pulse)
        seq.add_block(post_sl_pulse)
        if n < defs["n_pulses"] - 1:
            seq.add_block(pp.make_delay(defs["td"] - sl_time_per_sat))

    if FLAG_POST_PREP_SPOIL:
        seq.add_block(gx_spoil, gy_spoil, gz_spoil)
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
