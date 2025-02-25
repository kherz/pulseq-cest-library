from pathlib import Path
import numpy as np
import pypulseq as pp
from bmctool.utils.seq.write import write_seq

def generate_WASABI_sequence(time_min=11.5, time_max=15.5):
    # get id of generation file
    seqid = Path(__file__).stem + "_python"
    folder = Path.cwd()

    # general settings
    AUTHOR = "Patrick Schuenke"
    FLAG_PLOT_SEQUENCE = True  # plot preparation block?
    FLAG_CHECK_TIMING = False  # perform a timing check at the end of the sequence?
    FLAG_POST_PREP_SPOIL = True  # add spoiler after preparation block?

    # sequence definitions
    defs = {
        "b1pa": 3.7,  # B1 peak amplitude [µT]
        "b1rms": 3.7,  # B1 RMS amplitude [µT]
        "b0": 3,  # B0 [T]
        "n_pulses": 1,  # number of pulses  
        "tp": 5e-3,  # pulse duration [s]
        "trec": 3,  # recovery time [s]
        "trec_m0": 12,  # recovery time before M0 [s]
        "m0_offset": -300,  # m0 offset [ppm]
        "offsets_ppm": np.append(-300, np.linspace(-2, 2, 31)),  # offset vector [ppm]
        "spoiling": "1" if FLAG_POST_PREP_SPOIL else "0"
    }

    defs["num_meas"] = defs["offsets_ppm"].size
    defs["tsat"] = defs["tp"]
    defs["seq_id_string"] = seqid
    seq_filename = seqid + ".seq"

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
    defs["freq"] = defs["b0"] * GAMMA_HZ  # Larmor frequency [Hz]

    # Spoiler
    spoil_amp = 0.8 * sys.max_grad
    rise_time = 1.0e-3
    spoil_dur = 6.5e-3
    gx_spoil, gy_spoil, gz_spoil = [
        pp.make_trapezoid(channel=c, system=sys, amplitude=spoil_amp, duration=spoil_dur, rise_time=rise_time)
        for c in ["x", "y", "z"]
    ]

    # RF pulses
    flip_angle_sat = defs["b1pa"] * GAMMA_HZ * 2 * np.pi * defs["tp"]
    sat_pulse = pp.make_block_pulse(flip_angle=flip_angle_sat, duration=defs["tp"], system=sys)

    # Pseudo ADC event
    pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)

    # Delays
    trec_delay = pp.make_delay(defs["trec"])
    m0_delay = pp.make_delay(defs["trec_m0"])

    # Sequence object
    seq = pp.Sequence()

    # Offset calculations
    offsets_hz = defs["offsets_ppm"] * defs["freq"]

    for m, offset in enumerate(offsets_hz):
        print(f"#{m + 1} / {len(offsets_hz)} : offset {offset / defs['freq']:.2f} ppm ({offset:.3f} Hz)")

        if offset == defs["m0_offset"] * defs["freq"]:
            if defs["trec_m0"] > 0:
                seq.add_block(m0_delay)
        else:
            if defs["trec"] > 0:
                seq.add_block(trec_delay)

        sat_pulse.freq_offset = offset
        seq.add_block(sat_pulse)

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

    if FLAG_PLOT_SEQUENCE:
        seq.plot(time_range=(int(time_min), int(time_max)))

if __name__ == "__main__":
    generate_WASABI_sequence()
