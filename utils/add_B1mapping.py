# Add B1 Map Methode
# The sequence of saturated and unsaturated readouts defined here is the sequence implemented in ICE.
#
# doesn't support Pypulseq 1.4.5/1.4.6
# only for Pulseq >1.5.0!
#
# Martin Freudensprung 2026
# martin.freudensprung@fau.de


import numpy as np
import pypulseq as pp

system = pp.Opts(
    max_grad=40,
    grad_unit="mT/m",
    max_slew=130,
    slew_unit="T/m/s",
    rf_ringdown_time=30e-6,
    rf_dead_time=100e-6,
    rf_raster_time=1e-6,
    gamma=42576400,
)

def add_2scan_B1_Maps(seq, defs, system = system, shim_array = pp.get_tx_mode('Nova_Head_8Tx_CP')):


    tp_70deg = 2e-3  # s
    fa_70deg = np.pi * 70.0 / 180.0  # 70° to radians

    rfPulse_70deg = pp.make_block_pulse(flip_angle=fa_70deg,
                                           duration=tp_70deg,
                                           system=system)


    gx_spoil, gy_spoil, gz_spoil = [
        pp.make_trapezoid(channel=c,
                          system=system,
                          amplitude=0.8 * system.max_grad,  # Hz/m
                          duration=6.5e-3, # complete spoiler duration in seconds
                          rise_time=1.0e-3) # spoiler rise time in seconds
        for c in ["x", "y", "z"]
    ]

    if "trec_B1_map" not in defs:
        defs["trec_B1_map"] = 4

    rec_d = pp.make_delay(defs["trec_B1_map"])
    # pseudo ADC event
    pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

    #Add Pause before
    seq.add_block(rec_d)  # You never know the protocol history -> relax

    print("Add unsaturated image, excitation flip angle 0° (Mz,0)")
    # 0. Add unsaturated image, excitation flip angle 0° (Mz,0)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)

    seq.add_block(rec_d)
    # Add saturaed B1map Pulse
    print("Add saturated image, excitation flip angle 70° (Mz,70)")
    # 1. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg, pp.make_rf_shim(shim_array))
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    defs["numB1scans"] = 2
    defs["numB1Maps"] = 1
    defs["enableB1cor"] =1

    return seq, defs

def add_2scan_B1_Maps_rot(seq, defs, system = system, shim_array = pp.get_tx_mode('Nova_Head_8Tx_CP')):


    tp_70deg = 2e-3  # s
    fa_70deg = np.pi * 70.0 / 180.0  # 70° to radians

    rfPulse_70deg = pp.make_block_pulse(flip_angle=fa_70deg,
                                           duration=tp_70deg,
                                           system=system)


    gx_spoil, gy_spoil, gz_spoil = [
        pp.make_trapezoid(channel=c,
                          system=system,
                          amplitude=0.8 * system.max_grad,  # Hz/m
                          duration=6.5e-3, # complete spoiler duration in seconds
                          rise_time=1.0e-3) # spoiler rise time in seconds
        for c in ["x", "y", "z"]
    ]

    if "trec_B1_map" not in defs:
        defs["trec_B1_map"] = 4

    rec_d = pp.make_delay(defs["trec_B1_map"])
    # pseudo ADC event
    pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

    #Add Pause before
    seq.add_block(rec_d)  # You never know the protocol history -> relax

    print("Add saturated image, excitation flip angle 70° (Mz,70)")
    # 1. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg, pp.make_rf_shim(shim_array))
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    # Add saturaed B1map Pulse
    print("Add unsaturated image, excitation flip angle 0° (Mz,0)")
    # 0. Add unsaturated image, excitation flip angle 0° (Mz,0)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    defs["numB1scans"] = 2
    defs["numB1Maps"] = 1
    defs["enableB1cor"] =1

    return seq, defs

def add_4scan_B1_Maps(seq,
                    defs,
                    system = system,
                    shim_array_1 = pp.get_tx_mode('Nova_Head_8Tx_CP'),
                    shim_array_2 = pp.get_tx_mode('Nova_Head_8Tx_EP')):


    tp_70deg = 2e-3  # s
    fa_70deg = np.pi * 70.0 / 180.0  # 70° to radians

    rfPulse_70deg = pp.make_block_pulse(flip_angle=fa_70deg,
                                           duration=tp_70deg,
                                           system=system)


    gx_spoil, gy_spoil, gz_spoil = [
        pp.make_trapezoid(channel=c,
                          system=system,
                          amplitude=0.8 * system.max_grad,  # Hz/m
                          duration=6.5e-3, # complete spoiler duration in seconds
                          rise_time=1.0e-3) # spoiler rise time in seconds
        for c in ["x", "y", "z"]
    ]

    if "trec_B1_map" not in defs:
        defs["trec_B1_map"] = 4

    rec_d = pp.make_delay(defs["trec_B1_map"])
    # pseudo ADC event
    pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

    # Add Pause before
    seq.add_block(rec_d)  # You never know the protocol history -> relax

    print("Add unsaturated CP image, excitation flip angle 0° (Mz,0)")
    # 0. Add unsaturated image, excitation flip angle 0° (Mz,0)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    # first B1map
    print("Add saturated image, excitation flip angle 70° (Mz,70)")
    # 1. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg, pp.make_rf_shim(shim_array_1))
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    print("Add unsaturated image")
    # 0. Add unsaturated image, excitation flip angle 0° (Mz,0)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    # second B1map
    print("Add saturated image, excitation flip angle 70° (Mz,70)")
    # 1. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg, pp.make_rf_shim(shim_array_2))
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)


    defs["numB1scans"] = 4
    defs["numB1Maps"] = 2
    defs["enableB1cor"] =1

    return seq, defs

def add_3scan_B1_Maps(seq,
                    defs,
                    system = system,
                    shim_array_1 = pp.get_tx_mode('Nova_Head_8Tx_CP'),
                    shim_array_2 = pp.get_tx_mode('Nova_Head_8Tx_EP')):


    tp_70deg = 2e-3  # s
    fa_70deg = np.pi * 70.0 / 180.0  # 70° to radians

    rfPulse_70deg = pp.make_block_pulse(flip_angle=fa_70deg,
                                           duration=tp_70deg,
                                           system=system)

    gx_spoil, gy_spoil, gz_spoil = [
        pp.make_trapezoid(channel=c,
                          system=system,
                          amplitude=0.8 * system.max_grad,  # Hz/m
                          duration=6.5e-3, # complete spoiler duration in seconds
                          rise_time=1.0e-3) # spoiler rise time in seconds
        for c in ["x", "y", "z"]
    ]

    if "trec_B1_map" not in defs:
        defs["trec_B1_map"] = 4

    rec_d = pp.make_delay(defs["trec_B1_map"])
    # pseudo ADC event
    pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

    # Add Pause before
    seq.add_block(rec_d)  # You never know the protocol history -> relax

    print("Add unsaturated CP image, excitation flip angle 0° (Mz,0)")
    # 0. Add unsaturated image, excitation flip angle 0° (Mz,0)
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    # first B1map
    print("Add saturated image, excitation flip angle 70° (Mz,70)")
    # 1. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg, pp.make_rf_shim(shim_array_1))
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)
    seq.add_block(rec_d)

    # second B1map
    print("Add saturated image, excitation flip angle 70° (Mz,70)")
    # 1. Add unsaturated image, excitation flip angle 70° (Mz,70)
    seq.add_block(rfPulse_70deg, pp.make_rf_shim(shim_array_2))
    seq.add_block(gx_spoil, gy_spoil, gz_spoil)
    seq.add_block(pseudo_adc)

    defs["numB1scans"] = 3
    defs["numB1Maps"] = 2
    defs["enableB1cor"] = 1

    return seq, defs
