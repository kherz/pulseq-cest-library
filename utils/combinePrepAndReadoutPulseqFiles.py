# Skript to combine Pulseq-CEST Files with a Readout
# exchange the single ADC Event by the provided Readout

# To use the online Reconstruction at the scanner, the lables had to be set correctly
# This works only stable in Pulseq Version >1.5
# take care to only use Pulseq 1.5 implementations. Mixing different Version for CEST Preperation and Readout
# lead to problems

# June 2026
# Martin Freudensprung (Uniklinik Erlangen, Germany)
# martin.freudensprung@uk-erlangen.de

#%% ===========================================================================
# Imports
# =============================================================================
import copy
import pypulseq as pp

from matplotlib import pyplot as plt

system = pp.Opts(max_grad=80,
                    grad_unit='mT/m',
                    max_slew=180,
                    slew_unit='T/m/s',
                    rf_ringdown_time=20e-6,
                    rf_dead_time=100e-6,
                    adc_dead_time=10e-6)

#%% ===========================================================================
# Read CEST prep
# =============================================================================
seq_CEST = pp.Sequence()

seq_CEST.read(r"...")

#%% ===========================================================================
# Read Readout
seq_SNAPSHOT = pp.Sequence()
seq_SNAPSHOT.read(r'...\seq\2D_WIP_2mm_REP_rot.seq')

#%% ===========================================================================
# Combine seq files
# =============================================================================
final_seq_name = r"...combined.seq"
combinedSequence = pp.Sequence()
label_rep = 0
# %%

prep_list = []
for prepIdx in range(len(seq_CEST.block_events)):
    cBlock = copy.deepcopy(seq_CEST.get_block(prepIdx + 1))
    if cBlock.adc is None:
        combinedSequence.add_block(cBlock)
    else:
        for roIdx in range(len(seq_SNAPSHOT.block_events)):
            sBlock = copy.deepcopy(seq_SNAPSHOT.get_block(roIdx + 1))

            label_dict = getattr(sBlock, "label", None)

            if label_dict is not None:
                for entry in label_dict.values():
                    if getattr(entry, "label", None) == "REP":
                        value = getattr(entry, "value", None)
                        setattr(entry, "value", label_rep)

            combinedSequence.add_block(sBlock)

        label_rep += 1

for key, value in seq_CEST.definitions.items():
    combinedSequence.set_definition(key, value)
for key, value in seq_SNAPSHOT.definitions.items():
    combinedSequence.set_definition(key, value)

# %% Plot
combinedSequence.plot(show_rf_shim=False, )

#%%

combinedSequence.write(final_seq_name)
print("safed File: "+ str(final_seq_name))

