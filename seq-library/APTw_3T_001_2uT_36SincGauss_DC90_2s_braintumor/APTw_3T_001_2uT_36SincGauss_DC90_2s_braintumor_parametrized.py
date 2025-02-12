from pathlib import Path
import numpy as np
import pypulseq as pp
from bmctool.utils.pulses.calc_power_equivalents import calc_power_equivalent
from bmctool.utils.seq.write import write_seq

def generate_sequence(time_min=3, time_max=6):
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
        "b1pa": 1.78,  # B1 peak amplitude [µT]
        "b0": 3,  # B0 [T]
        "n_pulses": 36,  # number of pulses  
        "tp": 50e-3,  # pulse duration [s]
        "td": 5e-3,  # interpulse delay [s]
        "trec": 3.5,  # recovery time [s]
        "trec_m0": 3.5,  # recovery time before M0 [s]
        "m0_offset": -300,  # m0 offset [ppm]
        "offsets_ppm": np.append(-300, np.linspace(-4, 4, 33)),  # offset vector [ppm]
        "spoiling": "1" if FLAG_POST_PREP_SPOIL else "0"
    }
    
    defs["dcsat"] = defs["tp"] / (defs["tp"] + defs["td"])
    defs["num_meas"] = defs["offsets_ppm"].size
    defs["tsat"] = defs["n_pulses"] * (defs["tp"] + defs["td"]) - defs["td"]
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
    sat_pulse = pp.make_sinc_pulse(
        flip_angle=flip_angle_sat, duration=defs["tp"], system=sys, time_bw_product=2, apodization=0.15
    )
    
    defs["b1rms"] = calc_power_equivalent(rf_pulse=sat_pulse, tp=defs["tp"], td=defs["td"], gamma_hz=GAMMA_HZ)
    
    # Pseudo ADC event
    pseudo_adc = pp.make_adc(num_samples=1, duration=1e-3)
    
    # Delays
    post_spoil_delay = pp.make_delay(50e-6)
    td_delay = pp.make_delay(defs["td"])
    trec_delay = pp.make_delay(defs["trec"])
    m0_delay = pp.make_delay(defs["trec_m0"])
    
    # Sequence object
    seq = pp.Sequence()
    
    # Offset calculations
    offsets_hz = defs["offsets_ppm"] * defs["freq"]
    
    for m, offset in enumerate(offsets_hz):
        print(f"#{m + 1} / {len(offsets_hz)} : offset {offset / defs['freq']:.2f} ppm ({offset:.3f} Hz)")
        accum_phase = 0
        
        if offset == defs["m0_offset"] * defs["freq"]:
            if defs["trec_m0"] > 0:
                seq.add_block(m0_delay)
        else:
            if defs["trec"] > 0:
                seq.add_block(trec_delay)
        
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
    
    if FLAG_PLOT_SEQUENCE:
        seq.plot(time_range=(int(time_min), int(time_max)))

if __name__ == "__main__":
    generate_sequence()
