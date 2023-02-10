"""
Script to output a seq file for a T1 prep with an adiabatic hyperbolic secant preparation pulse.
"""

import os

import numpy as np
from bmctool.utils.pulses.make_hypsec_half_passage import make_hypsec_half_passage_rf
from bmctool.utils.seq.write import write_seq
from pypulseq import Opts
from pypulseq import Sequence
from pypulseq import make_adc
from pypulseq import make_delay
from pypulseq import make_trapezoid

# get id of generation file
seqid = os.path.splitext(os.path.basename(__file__))[0]

# general settings
author = 'Patrick Schuenke'
plot_sequence = False  # plot preparation block?
convert_to_1_3 = True  # convert seq-file to a pseudo version 1.3 file?
check_timing = True  # Perform a timing check at the end of the sequence

# inversion times (number of inversion times defines the number of measurements)
TI = np.array([10, 6, 5, 4, 3, 2.5, 2, 1.5, 1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1])

# sequence definitions (everything in seq_defs will be written to definitions of the .seq-file)
seq_defs: dict = {}
seq_defs['b1rms'] = 20
seq_defs['n_pulses'] = 3  # number of pulses  #
seq_defs['tp'] = 8e-3  # pulse duration [s]
seq_defs['trec'] = 1  # recovery time [s]
seq_defs['TI'] = TI  # inversion time vector
seq_defs['offsets_ppm'] = np.zeros(TI.shape)  # offset vector [ppm]

seq_defs['num_meas'] = seq_defs['offsets_ppm'].size  # number of repetition
seq_defs['seq_id_string'] = seqid  # unique seq id
seq_filename = seq_defs['seq_id_string'] + '.seq'

# scanner limits
sys = Opts(max_grad=40, grad_unit='mT/m', max_slew=130, slew_unit='T/m/s',
           rf_ringdown_time=30e-6, rf_dead_time=100e-6, rf_raster_time=1e-6)

gamma_hz = 42.5764

# ===========
# PREPARATION
# ===========

# spoilers
spoil_amp0 = 0.8 * sys.max_grad  # Hz/m
spoil_amp1 = -0.7 * sys.max_grad  # Hz/m
spoil_amp2 = 0.6 * sys.max_grad  # Hz/m

rise_time = 1.0e-3  # spoiler rise time in seconds
spoil_dur = 5.5e-3  # complete spoiler duration in seconds

gx_spoil0, gy_spoil0, gz_spoil0 = [make_trapezoid(channel=c, system=sys, amplitude=spoil_amp0, duration=spoil_dur,
                                                  rise_time=rise_time) for c in ['x', 'y', 'z']]
gx_spoil1, gy_spoil1, gz_spoil1 = [make_trapezoid(channel=c, system=sys, amplitude=spoil_amp1, duration=spoil_dur,
                                                  rise_time=rise_time) for c in ['x', 'y', 'z']]
gx_spoil2, gy_spoil2, gz_spoil2 = [make_trapezoid(channel=c, system=sys, amplitude=spoil_amp2, duration=spoil_dur,
                                                  rise_time=rise_time) for c in ['x', 'y', 'z']]

# RF pulses
hs_pulse = make_hypsec_half_passage_rf(amp=seq_defs['b1rms'], system=sys)

# ADC events
pseudo_adc = make_adc(num_samples=1, duration=1e-3)  # (not played out; just used to split measurements)

# DELAYS
trec_delay = make_delay(seq_defs['trec'])

# Sequence object
seq = Sequence()

# ===
# RUN
# ===

for m, t_prep in enumerate(TI):
    # print progress/offset
    print(f' {m + 1} / {len(TI)} : TI = {t_prep} s')

    # add delay
    seq.add_block(trec_delay)

    # add preparation (adiabatic excitation + spoiler) block
    for i in range(seq_defs['n_pulses']):
        seq.add_block(hs_pulse)
        if i % 3 == 0:
            seq.add_block(gx_spoil0, gy_spoil1, gz_spoil2)
        elif i % 2 == 0:
            seq.add_block(gx_spoil2, gy_spoil0, gz_spoil1)
        else:
            seq.add_block(gx_spoil1, gy_spoil2, gz_spoil0)

    # add variable inversion time delay
    seq.add_block(make_delay(t_prep))

    # add pseudo adc
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
