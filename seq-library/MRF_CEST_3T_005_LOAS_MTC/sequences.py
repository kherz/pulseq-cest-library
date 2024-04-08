import numpy as np
import pypulseq as pp
import os


def create_defs_loas_10():
    defs = {}
    defs['b1'] = [1.1, 1.9, 1.9, 0.9, 1.9, 1.1, 1.5, 1.9, 1.9, 0.9]
    defs['offsets_ppm'] = [9.1, 8.9, 50.0, 50.0, 11.5, 8.6, 20.7, 50.0, 50.0, 50.0]
    defs['tsat'] = [1.5, 2.0, 2.0, 0.4, 0.5, 2.0, 2.0, 2.0, 1.7, 0.4]
    defs['trec'] = [4.5, 3.5, 4.5, 3.5, 4.5, 3.9, 3.5, 4.4, 4.5, 3.5]
    return defs


def create_defs_loas_40():
    defs = {}
    defs['b1'] = [1.8, 1.2, 1.9, 1.9, 1.8, 1.3, 1.2, 1.9, 1.9, 0.9, 1.3, 1.4, 1.3, 1.2, 1.2, 1.1, 1.9, 0.9, 0.9, 1.5,
                        1.9, 1.9, 1.9, 1.9, 1.3, 1.1, 1.5, 1.2, 0.9, 1.9, 1.9, 0.9, 0.9, 1.8, 1.9, 1.4, 1.9, 1.2, 1.9, 1.9]
    defs['offsets_ppm'] = [50.0, 10.2, 23.6, 15.6, 11.3, 15.2, 8.9, 23.2, 28.4, 50.0, 28.6, 19.9, 11.9, 11.7, 9.7, 8.7, 9.0, 50.0,
                                49.9, 35.6, 50.0, 50.0, 10.1, 50.0, 10.1, 8.7, 35.1, 9.4, 49.8, 34.8, 10.4, 50.0, 50.0, 18.9, 8.1, 24.2,
                                10.8, 9.6, 50.0, 50.0]
    defs['tsat'] = [2.0, 2.0, 1.5, 0.7, 0.4, 2.0, 1.9, 2.0, 2.0, 0.4, 2.0, 1.9, 0.8, 2.0, 2.0, 1.9, 1.9, 0.4, 0.4, 2.0,
                        1.9, 1.9, 0.9, 0.4, 0.5, 2.0, 1.9, 2.0, 0.4, 0.6, 0.4, 0.4, 0.4, 0.7, 2.0, 2.0, 1.0, 1.0, 2.0, 2.0]
    defs['trec'] = [4.5, 3.7, 4.4, 4.5, 4.3, 3.6, 3.7, 3.5, 3.5, 3.5, 3.5, 3.6, 4.5, 3.8, 3.5, 4.2, 3.7, 3.5, 3.5, 4.5,
                        4.5, 4.5, 4.2, 3.5, 3.5, 3.9, 4.4, 3.6, 3.8, 4.4, 4.5, 3.5, 4.0, 4.5, 3.6, 3.5, 4.0, 4.4, 4.5, 4.2]

    return defs    


def create_defs_loas(N=10, gamma = 267.5153, freq = 127.7153, clinical=False):
    if N == 10:
        defs = create_defs_loas_10()
    elif N == 40:
        defs = create_defs_loas_40()
    else:
        raise ValueError(f"Invalid length: {N}")
    
    for key in ['b1', 'offsets_ppm', 'tsat', 'trec']:
        defs[key] = np.array(defs[key])

    if clinical:
        defs['tp'] = 100e-3
        defs['td'] = 100e-3
        defs['dcsat'] = defs['tp'] / (defs['tp'] + defs['td'])
        defs['n_pulses'] = np.ceil(defs['tsat'] / (defs['tp'] + defs['td'])).astype(int)

        lims = pp.Opts(
            max_grad=40,
            grad_unit="mT/m",
            max_slew=130,
            slew_unit="T/m/s",
            rf_ringdown_time=30e-6,
            rf_dead_time=100e-6,
            rf_raster_time=1e-6,
            gamma=gamma / 2 / np.pi * 1e6,
        )
        
        defs["gamma_hz"] = lims.gamma * 1e-6
        defs["freq"] = freq
        defs['B0'] = defs['freq'] / defs["gamma_hz"]

        defs['spoiling'] = True

        return defs, lims

    else:
        return defs


def write_sequence_clinical(defs: dict, seq_fn:str, lims = None):
    """
    Create clinical pulsed-wave sequence for CEST with complex readout
    :param defs: sequence definitions
    :param seq_fn: sequence filename
    :param lims: scanner limits
    :return: sequence object
    """

    GAMMA_HZ = defs["gamma_hz"]

    tp_sl = 1e-3  # duration of tipping pulse for sl
    td_sl = lims.rf_dead_time + lims.rf_ringdown_time  # delay between tip and sat pulse
    sl_time_per_sat = 2 * (tp_sl + td_sl)  # additional time of sl pulses for 1 sat pulse
    assert (defs["trec"] >= sl_time_per_sat).all(), "DC too high for SL preparation pulses!"

    sl_pause_time = 250e-6

    # spoiler
    spoil_amp = 0.8 * lims.max_grad  # Hz/m
    rise_time = 1.0e-3  # spoiler rise time in seconds
    spoil_dur = 4500e-6 + rise_time  # complete spoiler duration in seconds

    gx_spoil, gy_spoil, gz_spoil = [
        pp.make_trapezoid(channel=c, system=lims, amplitude=spoil_amp, duration=spoil_dur, rise_time=rise_time)
        for c in ["x", "y", "z"]
    ]
    sl_pause_time = sl_pause_time - lims.rf_dead_time - lims.rf_ringdown_time

    min_fa = 1

    pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)
    offsets_hz = defs["offsets_ppm"] * defs["freq"]  # convert from ppm to Hz

    phase_cycling = 50 / 180 * np.pi
    seq = pp.Sequence()

    for m, b1 in enumerate(defs["b1"]):
        # reset accumulated phase
        accum_phase = 0

        # prep and set rf pulse
        flip_angle_sat = b1 * GAMMA_HZ * 2 * np.pi * defs["tp"]
        sat_pulse = pp.make_block_pulse(flip_angle=flip_angle_sat, duration=defs["tp"], freq_offset=offsets_hz[m],
                                        system=lims)
        accum_phase = np.mod(offsets_hz[m] * 2 * np.pi * defs["tp"], 2 * np.pi)

        # prep spin lock pulses
        flip_angle_tip = np.arctan(b1 / (defs["offsets_ppm"][m] * defs["B0"] + 1e-8))
        pre_sl_pulse = pp.make_block_pulse(flip_angle=flip_angle_tip, duration=tp_sl, phase_offset=-(np.pi / 2),
                                           system=lims)
        post_sl_pulse = pp.make_block_pulse(flip_angle=flip_angle_tip, duration=tp_sl,
                                            phase_offset=accum_phase + (np.pi / 2), system=lims)

        sat_pulse.freq_offset = offsets_hz[m]
        for n in range(defs["n_pulses"][m]):
            pre_sl_pulse.phase_offset = pre_sl_pulse.phase_offset + phase_cycling
            sat_pulse.phase_offset = sat_pulse.phase_offset + phase_cycling
            post_sl_pulse.phase_offset = post_sl_pulse.phase_offset + phase_cycling

            if b1 == 0:
                seq.add_block(pp.make_delay(defs["tp"]))
                seq.add_block(pp.make_delay(sl_time_per_sat))
            else:
                if flip_angle_tip > min_fa/180*np.pi:
                    seq.add_block(pre_sl_pulse)
                    seq.add_block(pp.make_delay(sl_pause_time))
                else:
                    seq.add_block(pp.make_delay(pp.calc_duration(pre_sl_pulse)+sl_pause_time))

                seq.add_block(sat_pulse)

                if flip_angle_tip > min_fa/180*np.pi:
                    seq.add_block(pp.make_delay(sl_pause_time))
                    seq.add_block(post_sl_pulse)
                else:
                    seq.add_block(pp.make_delay(pp.calc_duration(post_sl_pulse)+sl_pause_time))

                if n < defs["n_pulses"][m] - 1:
                    seq.add_block(pp.make_delay(defs["td"] - sl_time_per_sat))

        seq.add_block(pp.make_delay(100e-6)) # hardware related delay

        if defs["spoiling"]:
            seq.add_block(gx_spoil, gy_spoil, gz_spoil)
            seq.add_block(pp.make_delay(100e-6)) # hardware related delay

        # Readout
        seq.add_block(pseudo_adc)
        readout_time = pp.calc_duration(pseudo_adc)

        # add delay
        if m < len(defs["b1"]) - 1:
            seq.add_block(pp.make_delay(defs["trec"][m]-readout_time))

    def_fields = defs.keys()
    for field in def_fields:
        seq.set_definition(field, defs[field])

    seq.write(seq_fn)
    return seq
