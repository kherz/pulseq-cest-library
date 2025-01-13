## utils

### 1. Pulseq-CEST merge to Pulseq
Pulseq-CEST files contain only a single ADC event after the saturation phase.
The script combinePrepAndReadoutPulseqFiles.m exchanges the single ADC event by the provided RO.
The script demoCombination.m gives and short exmaple.

#### Important information
- Total duration is removed from seq header since it could lead to problems during scanning
- Current pulseq version 1.4.2 has a bug in setBlock: Workaround as suggested here: https://github.com/pulseq/pulseq/issues/53
- It is possible to combined different pulseq versions


### 2. Arbitrary pulse-shape to Pulseq-CEST
A function / snippet to load arbitrary pulse shapes into a pulseq object.

