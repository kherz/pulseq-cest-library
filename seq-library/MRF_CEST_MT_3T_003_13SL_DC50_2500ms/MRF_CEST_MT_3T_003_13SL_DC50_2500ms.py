# MRF_CEST_MT_3T_003_13SL_DC50_2500ms
# Creates a sequence file for a MT CEST Fingerprinting protocol
# Reference: Perlman, O., Ito, H., Herz, K. et al. 
# Quantitative imaging of apoptosis following oncolytic virotherapy by magnetic resonance fingerprinting aided by deep learning. 
# Nat. Biomed. Eng, 2022, https://doi.org/10.1038/s41551-021-00809-7

from cmath import atan
import os
import numpy as np
from pypulseq.Sequence.sequence import Sequence
from pypulseq.make_adc import make_adc
from pypulseq.make_delay import make_delay
from pypulseq.make_trap_pulse import make_trapezoid
from pypulseq.make_block_pulse import make_block_pulse
from pypulseq.opts import Opts
from bmctool.utils.seq.write import write_seq

# get id of generation file
seqid = os.path.splitext(os.path.basename(__file__))[0]

# general settings
author = 'Kai Herz'
convert_to_1_3 = False  # convert seq-file to a version 1.3 file? Needed for pypulseq < v1.3.1 only!

# sequence definitions (everything in seq_defs will be written to definitions of the .seq-file)
num_meas = 30
seq_defs:dict = {}
seq_defs['num_meas'] = 30
seq_defs['b1cwpe']      = np.array([2.0, 2.0, 1.7, 1.5, 1.2, 1.2, 3.0, 0.5, 3.0, 1.0, 2.2, 3.2, 1.5, 0.7, 1.5, 2.2, 2.5, 1.2, 3.0, 0.2, 1.5, 2.5, 0.7, 4.0, 3.2, 3.5, 1.5, 2.7, 0.7, 0.5])
seq_defs['offsets_ppm'] = np.array([8,   6,   6,   10,  10,  10,  8,   6,   8,   14,  14,  10,  6,   10,  8,   6,   8,   10,  14,  14,  6,   8,   14,  6,   14,  14,  14,  8,   10,  8  ])
tr = np.ones(seq_defs['num_meas'])*3.5
seq_defs['tsat'] = np.ones(seq_defs['num_meas'])*2.5
seq_defs['trec'] = tr - seq_defs['tsat']
seq_defs['b0'] = 3  # B0 [T]
seq_defs['n_pulses'] = 13  # number of pulses  #
seq_defs['tp'] = 100e-3  # pulse duration [s]
seq_defs['td'] = 100e-3  # interpulse delay [s]
seq_defs['dcsat'] = (seq_defs['tp']) / (seq_defs['tp'] + seq_defs['td'])  # duty cycle
seq_defs['seq_id_string'] = seqid  # unique seq id
seq_filename = seq_defs['seq_id_string'] + '.seq'

# scanner limits
sys = Opts(max_grad=40, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
           rf_ringdown_time=50e-6, rf_dead_time=200e-6, rf_raster_time=1e-6)

gamma_hz = 42.5764

# ===========
# PREPARATION
# ===========

# spin lock specific preparation
# td should be between block pulses, so we fit the tipping pulses in the
# delay between pulses
tp_sl        = 1e-3                                  # duration of tipping pulse for sl
td_sl        = (sys.rf_dead_time + sys.rf_ringdown_time) # delay between tip and sat pulse
sl_time_per_sat = 2*(tp_sl+td_sl) # additional time of sl pulses for 1 sat pulse
if any(seq_defs['trec']  < sl_time_per_sat):
    print('DC too high for SL prepatration pulses!')
    quit()

# spoiler
spoil_amp = 0.8 * sys.max_grad  # Hz/m
rise_time = 1.0e-3  # spoiler rise time in seconds
spoil_dur = 6.5e-3  # complete spoiler duration in seconds

gx_spoil, gy_spoil, gz_spoil = [make_trapezoid(channel=c, system=sys, amplitude=spoil_amp, duration=spoil_dur,
                                               rise_time=rise_time) for c in ['x', 'y', 'z']]

# ADC events
pseudo_adc = make_adc(num_samples=1, duration=1e-3)  # (not played out; just used to split measurements)

# Sequence object
seq = Sequence()

# ===
# RUN
# ===

offsets_hz = seq_defs['offsets_ppm'] * gamma_hz * seq_defs['b0']  # convert from ppm to Hz

for m, b1 in enumerate(seq_defs['b1cwpe']):

    
    # add delay
    seq.add_block(make_delay(seq_defs['trec'][m]))

    # prep sat_pulse
    flip_angle_sat = b1 * gamma_hz * 2 * np.pi * seq_defs['tp']
    sat_pulse = make_block_pulse(flip_angle=flip_angle_sat, freq_offset=offsets_hz[m], duration=seq_defs['tp'], system=sys) #
    accum_phase = np.mod(offsets_hz[m]*2*np.pi*seq_defs['tp'],2*np.pi)
    # prep tipping pulses
    flip_angle_tip = atan(b1/(seq_defs['offsets_ppm'][m]*seq_defs['b0']))
    pre_sl_pulse  = make_block_pulse(flip_angle=flip_angle_tip, duration=tp_sl, phase_offset=-(np.pi/2), system=sys) #
    post_sl_pulse = make_block_pulse(flip_angle=flip_angle_tip, duration=tp_sl, phase_offset=accum_phase+(np.pi/2), system=sys) #
    # sl phase cycling
    phase_cycling = 50/180*np.pi

    for n in range(seq_defs['n_pulses']):
        pre_sl_pulse.phase_offset = pre_sl_pulse.phase_offset + phase_cycling
        sat_pulse.phase_offset = sat_pulse.phase_offset + phase_cycling
        post_sl_pulse.phase_offset = post_sl_pulse.phase_offset + phase_cycling
        seq.add_block(pre_sl_pulse)
        seq.add_block(sat_pulse)
        seq.add_block(post_sl_pulse)
        if n < seq_defs['n_pulses']-1:
            seq.add_block(make_delay(seq_defs['td']-sl_time_per_sat))

    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)

write_seq(seq=seq,
          seq_defs=seq_defs,
          filename=seq_filename,
          author=author,
          use_matlab_names=True,
          convert_to_1_3=convert_to_1_3)
