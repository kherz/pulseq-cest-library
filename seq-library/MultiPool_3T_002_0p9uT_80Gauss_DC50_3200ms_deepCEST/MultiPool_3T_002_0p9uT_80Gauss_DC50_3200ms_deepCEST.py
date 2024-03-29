# MultiPool_3T_002_0p9uT_80Gauss_DC50_3200ms_deepCEST
# Creates a sequence file for a multi pool protocol for 3T as used for the 3T to 7T deepCEST pipeline by Leonie Hunger
# this is the 0.9 µT protocol from MultiPool_3T_001_0p9uT_80Gauss_DC50_3200ms_GLINT
but with a slightly
# different offsetlist compared to the original publications.

import os
import numpy as np
from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.make_gauss_pulse import make_gauss_pulse
from pypulseq.opts import Opts
from bmctool.utils.pulses.calc_power_equivalents import calc_power_equivalent
from bmctool.utils.seq.write import write_seq

# get id of generation file
seqid = os.path.splitext(os.path.basename(__file__))[0]

# general settings
author = 'Lukas Kamm'
plot_sequence = True  # plot preparation block?
convert_to_1_3 = False  # convert seq-file to a version 1.3 file? Needed for pypulseq < v1.3.1 only!

# sequence definitions (everything in seq_defs will be written to definitions of the .seq-file)
b1 : float = 0.6 # B1 peak amplitude [µT] (the cw power equivalent will be calculated and written to seq_defs below)
seq_defs:dict = {}
seq_defs['b0'] = 3  # B0 [T]
seq_defs['n_pulses'] = 80  # number of pulses  #
seq_defs['tp'] = 20.48e-3  # pulse duration [s]
seq_defs['td'] = 20e-3  # interpulse delay [s]
seq_defs['trec'] = 0  # recovery time [s]
seq_defs['trec_m0'] = 12  # recovery time before M0 [s]
seq_defs['m0_offset'] = -300  # m0 offset [ppm]
seq_defs['offsets_ppm'] = np.append(
    seq_defs['m0_offset'], np.array([-100,  -50,   -40,   -30,   -20,   -10,    -9.5,   -9,    -8.5,   -8,    -7.5,   -7,    -6.5,   -6,    -5.5,   -5,    -4.5,   -4,    -3.5,   -3,    -2.5,   -2,    -1.5,   -1,    -0.6,   -0.4,   -0.2,    0,     0.2,    0.4,     0.6,    1,     1.5,    2,     2.5,    3,     3.5,    4,     4.5,    5,     5.5,    6,     6.5,    7,     7.5,    8,     8.5,    9,     9.5,   10,    20,    30,    40,    50,   100]))

seq_defs['dcsat'] = (seq_defs['tp']) / (seq_defs['tp'] + seq_defs['td'])  # duty cycle
seq_defs['num_meas'] = seq_defs['offsets_ppm'].size  # number of repetition
seq_defs['tsat'] = seq_defs['n_pulses'] * (seq_defs['tp'] + seq_defs['td']) - seq_defs['td']  # saturation time [s]
seq_defs['seq_id_string'] = seqid  # unique seq id

seq_filename = seq_defs['seq_id_string'] + '.seq'

# scanner limits
sys = Opts(max_grad=40, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
           rf_ringdown_time=30e-6, rf_dead_time=100e-6, rf_raster_time=1e-6)

gamma_hz = 42.5764

# ===========
# PREPARATION
# ===========

# spoiler
spoil_amp = 0.8 * sys.max_grad  # Hz/m
rise_time = 1.0e-3  # spoiler rise time in seconds
spoil_dur = 6.5e-3  # complete spoiler duration in seconds

gx_spoil, gy_spoil, gz_spoil = [make_trapezoid(channel=c, system=sys, amplitude=spoil_amp, duration=spoil_dur,
                                               rise_time=rise_time) for c in ['x', 'y', 'z']]

# RF pulses
flip_angle_sat = b1 * gamma_hz * 2 * np.pi * seq_defs['tp']
sat_pulse = make_gauss_pulse(flip_angle=flip_angle_sat, duration=seq_defs['tp'], system=sys,
                             time_bw_product=0.2, apodization=0.5)  # siemens-like gauss

seq_defs['b1cwpe'] = calc_power_equivalent(rf_pulse=sat_pulse, tp=seq_defs['tp'], td=seq_defs['td'], gamma_hz=gamma_hz)

# ADC events
pseudo_adc = make_adc(num_samples=1, duration=1e-3)  # (not played out; just used to split measurements)

# DELAYS
post_spoil_delay = make_delay(50e-6)
td_delay = make_delay(seq_defs['td'])
trec_delay = make_delay(seq_defs['trec'])
m0_delay = make_delay(seq_defs['trec_m0'])

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
    if offset == seq_defs['m0_offset'] * gamma_hz * seq_defs['b0']:
        if seq_defs['trec_m0'] > 0:
            seq.add_block(m0_delay)
    else:
        if seq_defs['trec'] > 0:
            seq.add_block(trec_delay)

    # set sat_pulse
    sat_pulse.freq_offset = offset
    for n in range(seq_defs['n_pulses']):
        sat_pulse.phase_offset = accum_phase % (2 * np.pi)
        seq.add_block(sat_pulse)
        accum_phase = (accum_phase + offset * 2 * np.pi * np.sum(np.abs(sat_pulse.signal) > 0) * 1e-6) % (2 * np.pi)
        if n < seq_defs['n_pulses']-1:
            seq.add_block(td_delay)

    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)

write_seq(seq=seq,
          seq_defs=seq_defs,
          filename=seqid+'.seq',
          author=author,
          use_matlab_names=True,
          convert_to_1_3=convert_to_1_3)

# plot the sequence
if plot_sequence:
    seq.plot(time_range=[seq_defs['trec_m0'], seq_defs['trec_m0']+seq_defs['tsat']])  # to plot all offsets, remove time_range argument
