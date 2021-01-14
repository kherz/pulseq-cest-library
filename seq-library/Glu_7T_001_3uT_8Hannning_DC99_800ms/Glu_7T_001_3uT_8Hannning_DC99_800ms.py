# Glu_7T_001_3uT_8Hannning_DC99_800ms_braintumor
# Creates a sequence file for an glutamate weighted protocol according to
# https://doi.org/10.1002/mrm.27362
#
# Kerstin Heinecke 2020
# kerstin.heinecke@ptb.de

import os
import numpy as np
from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.opts import Opts
from sim.utils.calc_power_equivalents import calc_power_equivalent
from sim.utils.seq.write_seq import write_seq
from sim.utils.seq.make_hanning import make_hanning_pulse

# get id of generation file
seqid = os.path.splitext(os.path.basename(__file__))[0]

# general settings
author = 'Kerstin Heinecke'
plot_sequence = False  # plot preparation block?
convert_to_1_3 = True  # convert seq-file to a pseudo version 1.3 file?

# sequence definitions (everything in seq_defs will be written to definitions of the .seq-file)
b1: float = 2.0  # B1 peak amplitude [µT] (the cw power equivalent will be calculated and written to seq_defs below)
seq_defs:dict = {}
seq_defs['b0'] = 7  # B0 [T]
seq_defs['n_pulses'] = 8  # number of pulses  #
seq_defs['tp'] = 99.8e-3  # pulse duration [s]
seq_defs['td'] = 0.2e-3  # interpulse delay [s]
seq_defs['trec'] = 15  # recovery time [s]
seq_defs['offsets_ppm'] = np.concatenate([[-100, -20], np.arange(-4.2, -1.8, 0.2),
                                          np.arange(1.8, 4.2, 0.2), [20, 100]])  # offset vector [ppm]
seq_defs['dcsat'] = seq_defs['tp'] / (seq_defs['tp'] + seq_defs['td'])  # duty cycle
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
flip_angle_sat = 3772 * np.pi / 180
sat_pulse = make_hanning_pulse(flip_angle_sat=flip_angle_sat, duration=seq_defs['tp'], system=sys)
seq_defs['b1cwpe'] = calc_power_equivalent(rf_pulse=sat_pulse, tp=seq_defs['tp'], td=seq_defs['td'], gamma_hz=gamma_hz)


# ADC events
pseudo_adc = make_adc(num_samples=1, duration=1e-3)  # (not played out; just used to split measurements)

# DELAYS
td_delay = make_delay(seq_defs['td'])
trec_delay = make_delay(seq_defs['trec'])

# Sequence object
seq = Sequence()
seq.shape_library

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
    seq.plot(time_range=[0, seq_defs['trec_m0']+seq_defs['tsat']])  # to plot all offsets, remove time_range argument
