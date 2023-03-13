# MultiPool_9p4T_001_multiB1sat_34Gauss_DC50_3p4s_brain
# Multi pool protocol for 9.4T according to
# https://cest-sources.org/doku.php?id=standard_cest_protocols
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
from bmctool.utils.pulses.calc_power_equivalents import calc_power_equivalent
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
defs["b1pa"] = np.array([0.6, 0.9, 1.1, 1.4])  # B1 peak amplitude [µT] (b1rms calculated below)
defs["freq"] = 400.2182
defs["n_pulses"] = 34  # number of pulses  #
defs["tp"] = 50e-3  # pulse duration [s]
defs["td"] = 50e-3  # interpulse delay [s]
defs["trec"] = 3.5  # recovery time [s]
defs["trec_m0"] = 12  # recovery time before M0 [s]
defs["m0_offset"] = -300  # m0 offset [ppm]

offset_list = [
    -50,
    -30,
    -20,
    -12,
    -8,
    -6,
    -5.6,
    -5.3,
    -5,
    -4.7,
    -4.4,
    -4.1,
    -3.8,
    -3.5,
    -3.2,
    -2.9,
    -2.6,
    -2.3,
    -2.1,
    -2,
    -1.7,
    -1.4,
    -1.15,
    -1,
    -0.85,
    -0.75,
    -0.5,
    -0.2,
    0,
    0.2,
    0.5,
    0.6,
    0.7,
    0.8,
    1.0,
    1.15,
    1.3,
    1.45,
    1.6,
    1.75,
    1.9,
    2.0,
    2.1,
    2.25,
    2.4,
    2.55,
    2.7,
    2.85,
    3.0,
    3.15,
    3.3,
    3.4,
    3.5,
    3.6,
    3.7,
    3.85,
    4.0,
    4.15,
    4.3,
    4.45,
    4.6,
    4.75,
    4.9,
    5.0,
    5.15,
    5.3,
    5.45,
    6.25,
    8.0,
    12.0,
    20.0,
    30.0,
    50.0,
]
defs["offsets_ppm"] = np.append(defs["m0_offset"], np.array(offset_list))

defs["dcsat"] = (defs["tp"]) / (defs["tp"] + defs["td"])  # duty cycle
defs["num_meas"] = defs["offsets_ppm"].size  # number of repetition
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
    rf_ringdown_time=30e-6,
    rf_dead_time=100e-6,
    rf_raster_time=1e-6,
    gamma=42576400,
)

GAMMA_HZ = sys.gamma * 1e-6
defs["b0"] = defs["freq"] / GAMMA_HZ  # magnetic field strength [T]

# ===========
# PREPARATION
# ===========

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
td_delay = pp.make_delay(defs["td"])
trec_delay = pp.make_delay(defs["trec"])
m0_delay = pp.make_delay(defs["trec_m0"])

# Sequence object
seq = pp.Sequence()

# ===
# RUN
# ===

defs["b1rms"] = []
offsets_hz = defs["offsets_ppm"] * defs["freq"]  # convert from ppm to Hz

for b1pa in defs["b1pa"]:
    # RF pulses
    flip_angle_sat = b1pa * GAMMA_HZ * 2 * np.pi * defs["tp"]
    sat_pulse = pp.make_gauss_pulse(
        flip_angle=flip_angle_sat, duration=defs["tp"], system=sys, time_bw_product=0.2, apodization=0.5
    )

    b1rms = calc_power_equivalent(rf_pulse=sat_pulse, tp=defs["tp"], td=defs["td"], gamma_hz=GAMMA_HZ)
    defs["b1rms"].append(b1rms)

    for m, offset in enumerate(offsets_hz):
        # print progress/offset
        print(f"#{m + 1} / {len(offsets_hz)} : offset {offset / defs['freq']:.2f} ppm ({offset:.3f} Hz)")

        # reset accumulated phase
        accum_phase = 0

        # add delay
        if offset == defs["m0_offset"] * defs["freq"]:
            if defs["trec_m0"] > 0:
                seq.add_block(m0_delay)
        else:
            if defs["trec"] > 0:
                # seq.add_block(trec_delay)
                pass

        # set sat_pulse
        sat_pulse.freq_offset = offset
        for n in range(defs["n_pulses"]):
            sat_pulse.phase_offset = accum_phase % (2 * np.pi)
            seq.add_block(sat_pulse)
            accum_phase = (accum_phase + offset * 2 * np.pi * np.sum(np.abs(sat_pulse.signal) > 0) * 1e-6) % (2 * np.pi)
            if n < defs["n_pulses"] - 1:
                seq.add_block(td_delay)
        
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
