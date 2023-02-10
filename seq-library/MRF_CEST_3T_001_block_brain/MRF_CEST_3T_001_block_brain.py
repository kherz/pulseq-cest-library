# MRF_CEST_3T_001_block_brain
# Creates a sequence file for an MRF-CEST protocol according to Figure 1 of:
# doi:10.1016/j.neuroimage.2019.01.034.
#
# Patrick Schuenke 2022
# patrick.schuenke@ptb.de

import os

import numpy as np
from bmctool.utils.seq.write import write_seq
from pypulseq import Opts
from pypulseq import Sequence
from pypulseq import make_adc
from pypulseq import make_block_pulse
from pypulseq import make_delay
from pypulseq import make_trapezoid

# get id of generation file
seqid = os.path.splitext(os.path.basename(__file__))[0]

# general settings
author = 'Patrick Schuenke'
plot_sequence = False  # plot preparation block?
convert_to_1_3 = False  # convert seq-file to a version 1.3 file? Needed for pypulseq < v1.3.1 only!
check_timing = True  # Perform a timing check at the end of the sequence

# sequence definitions (everything in seq_defs will be written to definitions of the .seq-file)
seq_defs: dict = {}
seq_defs['B1pa'] = np.array([1.2, 0.8, 2, 3])  # B1 peak amplitude [µT]
seq_defs['b0'] = 3  # B0 [T]
seq_defs['n_pulses'] = np.array([4, 2, 4, 3])  # number of pulses  #
seq_defs['tp'] = 200e-3  # pulse duration [s]
seq_defs['td'] = 10e-3  # interpulse delay [s]
seq_defs['trec'] = 2.5  # recovery time [s]
seq_defs['offsets_ppm'] = np.array([4, 3, 3.5, 10])  # offset vector [ppm]

seq_defs['dcsat'] = (seq_defs['tp']) / (seq_defs['tp'] + seq_defs['td'])  # duty cycle
seq_defs['num_meas'] = seq_defs['offsets_ppm'].size  # number of repetition
seq_defs['tsat'] = seq_defs['n_pulses'] * (seq_defs['tp'] + seq_defs['td']) - seq_defs['td']  # saturation time [s]
seq_defs['seq_id_string'] = seqid  # unique seq id

seq_filename = seq_defs['seq_id_string'] + '.seq'

# scanner limits
sys = Opts(max_grad=40, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
           rf_ringdown_time=30e-6, rf_dead_time=100e-6, rf_raster_time=1e-6)

gamma_hz =seq.sys.gamma*10e-6

# ===========
# PREPARATION
# ===========

# spoiler
spoil_amp = 0.8 * sys.max_grad  # Hz/m
rise_time = 1.0e-3  # spoiler rise time in seconds
spoil_dur = 6.5e-3  # complete spoiler duration in seconds

gx_spoil, gy_spoil, gz_spoil = [make_trapezoid(channel=c, system=sys, amplitude=spoil_amp, duration=spoil_dur,
                                               rise_time=rise_time) for c in ['x', 'y', 'z']]

# ADC events
pseudo_adc = make_adc(num_samples=1, duration=1e-3)  # (not played out; just used to split measurements)

# DELAYS
post_spoil_delay = make_delay(50e-6)
td_delay = make_delay(seq_defs['td'])
trec_delay = make_delay(seq_defs['trec'])

# Sequence object
seq = Sequence()

# ===
# RUN
# ===

offsets_hz = seq_defs['offsets_ppm'] * gamma_hz * seq_defs['b0']  # convert from ppm to Hz

for m, offset in enumerate(offsets_hz):
    # print progress/offset
    print(f' {m + 1} / {len(offsets_hz)} : offset {offset}')

    # reset accumulated phase
    accum_phase = 0

    # add delay
    seq.add_block(trec_delay)

    # prep and set rf pulse
    flip_angle_sat = seq_defs['B1pa'][m] * gamma_hz * 2 * np.pi * seq_defs['tp']
    sat_pulse = make_block_pulse(flip_angle=flip_angle_sat, duration=seq_defs['tp'], system=sys)
    sat_pulse.freq_offset = offset

    for n in range(seq_defs['n_pulses'][m]):
        sat_pulse.phase_offset = accum_phase % (2 * np.pi)
        seq.add_block(sat_pulse)
        accum_phase = (accum_phase + offset * 2 * np.pi * np.sum(np.abs(sat_pulse.signal) > 0) * 1e-6) % (2 * np.pi)
        if n < seq_defs['n_pulses'][m] - 1:
            seq.add_block(td_delay)

    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)

if check_timing:
    ok, error_report = seq.check_timing()
    if ok:
        print('\nTiming check passed successfully')
    else:
        print('\nTiming check failed! Error listing follows\n')
        print(error_report)

write_seq(seq=seq,
          seq_defs=seq_defs,
          filename=seq_filename,
          author=author,
          use_matlab_names=True,
          convert_to_1_3=convert_to_1_3)

# plot the sequence
if plot_sequence:
    seq.plot()  # to plot all offsets, remove time_range argument
