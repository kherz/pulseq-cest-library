# MultiPool_7T_001_0p57uT_120Gauss_DC60_3s_brain
# Multi pool protocol for 7T according to
# https://cest-sources.org/doku.php?id=standard_cest_protocols
#
#
# Tested with pypulseq version 1.4.6 and bmctool version 0.6.0
#
# created by
# Patrick Schuenke 2023
# adapted by
# Martin Freudensprung
# martin.freudesprung@fau.de

from pathlib import Path

import numpy as np
import pypulseq as pp
from bmctool.utils.pulses.calc_power_equivalents import calc_power_equivalent
from bmctool.utils.seq.write import write_seq
import os


addB1_Maps = True

# get id of generation file
seqid = Path(__file__).stem + "_python"

# get folder of generation file
folder = Path(__file__).parent

# general settings
AUTHOR = "Martin Freudensprung"
FLAG_PLOT_SEQUENCE = True  # plot preparation block?
FLAG_CHECK_TIMING = False  # perform a timing check at the end of the sequence?
FLAG_POST_PREP_SPOIL = True  # add spoiler after preparation block?

# sequence definitions
defs: dict = {}
defs["pulse_mode"]= 'MIMOSA'
defs["b1pa"] = 1.00  # B1 peak amplitude [µT] (b1rms calculated below)
defs["b0"] = 7  # B0 [T]
defs["n_pulses"] = 120  # number of pulses  #
defs["tp"] = 15e-3  # pulse duration [s]
defs["td"] = 10e-3  # interpulse delay [s]
defs["trec"] = 1  # recovery time [s]
defs["trec_m0"] = 12  # recovery time before M0 [s]
defs["m0_offset"] = -300  # m0 offset [ppm]
defs["n_dummies"] = 10
defs["offsets_ppm"] =  np.array(
            [
                -201, #1
                -201, #2
                -201, #3
                -201, #4
                -201, #5
                -201, #6
                -201, #7
                -201, #8
                -201, #9
                -201, #10
                -100,
                -50,
                -20,
                -12,
                -9,
                -7.25,
                -6.25,
                -5.5,
                -4.7,
                -4,
                -3.3,
                -2.7,
                -2,
                -1.7,
                -1.5,
                -1.1,
                -0.9,
                -0.6,
                -0.4,
                0,
                0.4,
                0.6,
                0.95,
                1.1,
                1.25,
                1.4,
                1.55,
                1.7,
                1.85,
                2,
                2.15,
                2.3,
                2.45,
                2.6,
                2.75,
                2.9,
                3.05,
                3.2,
                3.35,
                3.5,
                3.65,
                3.8,
                3.95,
                4.1,
                4.25,
                4.4,
                4.7,
                5.25,
                6.25,
                8,
                12,
                20,
                50,
                100,
            ])




defs["dcsat"] = (defs["tp"]) / (defs["tp"] + defs["td"])  # duty cycle
defs["num_meas"] = defs["offsets_ppm"].size  # number of repetition
if addB1_Maps:
    defs["num_meas"]+=4 # extra for for B1 mapping EP CP
defs["tsat"] = defs["n_pulses"] * (defs["tp"] + defs["td"]) - defs["td"]  # saturation time [s]
defs["seq_id_string"] = seqid  # unique seq id
defs["spoiling"] = "1" if FLAG_POST_PREP_SPOIL else "0"

defs['B1map_info'] = 'none_CP70deg_noneEP70deg'
defs["numB1scans"] = 2
defs["numB1Maps"] = 1
defs["enableB1cor"] = 1

seq_filename = defs["seq_id_string"] + ".seq"

# scanner limits Terra.X
sys = pp.Opts(max_grad=60,
        grad_unit='mT/m',
        max_slew=130,
        slew_unit='T/m/s',
        rf_ringdown_time=20e-6,
        rf_dead_time=100e-6,
        adc_dead_time=10e-6,
        B0=7.0,
        gamma=42576400)


GAMMA_HZ = sys.gamma * 1e-6
defs["freq"] = defs["b0"] * GAMMA_HZ  # Larmor frequency [Hz]

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

# RF pulses
flip_angle_sat = defs["b1pa"] * GAMMA_HZ * 2 * np.pi * defs["tp"]
sat_pulse = pp.make_gauss_pulse(
    flip_angle=flip_angle_sat, duration=defs["tp"], system=sys, time_bw_product=0.2, apodization=0.5)

if 1: #MIMOSA pulses
#cp mode
    sat_pulse_cp = pp.make_gauss_pulse(
        flip_angle=flip_angle_sat, duration=defs["tp"], system=sys, time_bw_product=0.2, apodization=0.5, shim_array=pp.set_tx_mode.CP_8ch())#ep mode
    sat_pulse_ep = pp.make_gauss_pulse(
        flip_angle=flip_angle_sat, duration=defs["tp"], system=sys, time_bw_product=0.2, apodization=0.5, shim_array=pp.set_tx_mode.EP_8ch())
else:  # fake CP
    sat_pulse_cp = pp.make_gauss_pulse(
        flip_angle=flip_angle_sat, duration=defs["tp"], system=sys, time_bw_product=0.2, apodization=0.5)
    sat_pulse_ep = pp.make_gauss_pulse(
        flip_angle=flip_angle_sat, duration=defs["tp"], system=sys, time_bw_product=0.2, apodization=0.5)

defs["b1rms"] = calc_power_equivalent(rf_pulse=sat_pulse, tp=defs["tp"], td=defs["td"], gamma_hz=GAMMA_HZ)
# pseudo ADC event
pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)
#pseudo_adc = pp.make_adc(num_samples=1, duration=5)

# delays
td_delay = pp.make_delay(defs["td"])
trec_delay = pp.make_delay(defs["trec"])
m0_delay = pp.make_delay(defs["trec_m0"])

# add B1 mapping for CP and EP with B1 presat method
tp_70deg = 2e-3  # s
fa_70deg = np.pi * 70.0 / 180.0  # 70° to radians
rfPulse_70deg = pp.make_block_pulse(flip_angle=fa_70deg, duration=tp_70deg, system=sys)
rec_d = pp.make_delay(4)


# Sequence object
seq = pp.Sequence()

# ===
# RUN
# ===

#CEST Offsets:

offsets_hz = defs["offsets_ppm"] * defs["freq"]  # convert from ppm to Hz
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
            seq.add_block(trec_delay)

    # set sat_pulse
    sat_pulse.freq_offset = offset
    sat_pulse_cp.freq_offset = offset
    sat_pulse_ep.freq_offset = offset
    for n in range(defs["n_pulses"]):
    #alternating pulses cp-ep-cp-ep-cp-ep-.....
        if n % 2 ==0: #even, cp mode
            sat_pulse_cp.phase_offset = accum_phase % (2 * np.pi)
            seq.add_block(sat_pulse_cp)
            accum_phase = (accum_phase + offset * 2 * np.pi * np.sum(np.abs(sat_pulse_cp.signal) > 0) * 1e-6) % (2 * np.pi)
        else:# odd, ep mode
            sat_pulse_ep.phase_offset = accum_phase % (2 * np.pi)
            seq.add_block(sat_pulse_ep)
            accum_phase = (accum_phase + offset * 2 * np.pi * np.sum(np.abs(sat_pulse_ep.signal) > 0) * 1e-6) % (2 * np.pi) 
        if n < defs["n_pulses"] - 1:
            seq.add_block(td_delay)
        print(accum_phase*(180/np.pi))
    if FLAG_POST_PREP_SPOIL:
        seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    
#B1 Maps
if addB1_Maps:
    # first B1map_CP (TODO)
    print("Add unsaturated CP image, excitation flip angle 0° (Mz,0)")
    # 0. Add unsaturated image, excitation flip angle 0° (Mz,0)
    seq.add_block(rec_d)  # You never know the protocol history -> relax
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    print("Add saturated CP image, excitation flip angle 70° (Mz,70)")
    # 1. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    # now switch to B1map_EP (TODO)
    print("Add unsaturated EP image, excitation flip angle 0° (Mz,0)")
    # 3. Add unsaturated image, excitation flip angle 0° (Mz,0)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    print("Add saturated EP image, excitation flip angle 70° (Mz,70)")
    # 4. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)


if FLAG_CHECK_TIMING:
    ok, error_report = seq.check_timing()
    if ok:
        print("\nTiming check passed successfully")
    else:
        print("\nTiming check failed! Error listing follows\n")
        print(error_report)

for key, value in defs.items():
    seq.set_definition(key, value)

#get the folder of the python file and save the seq file there
path = os.getcwd()
seq.write(path + "\\" + seqid + "_ptx" +'.seq')


seq.plot(time_range=(200,250))
#seq.plot()


