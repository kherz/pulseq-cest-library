function combinedSequence = combinePrepAndReadoutPulseqFiles(prep_fn, readout_fn, combined_fn)
%% combinePrepAndReadoutPulseqFiles
% script to exchange single ADC event by the proveded readout
% Inputs:
%   prep_fn - Path of CEST preparation with single ADC event.
%   readout_fn - Path of readout sequence
%   combined_fn - Path where the combined sequence should be saved
%   return
%   combinedSequence - Return of seq object of the combined pulseq seq
%   object
%
%% prepare sequence objects
prepSequence = mr.Sequence();
prepSequence.read(prep_fn);

readoutSequence = mr.Sequence();
readoutSequence.read(readout_fn);

combinedSequence = mr.Sequence();

%% loop through prep block
for prepIdx = 1:length(prepSequence.blockEvents)
    
    cBlock = prepSequence.getBlock(prepIdx);
    if isempty(cBlock.adc)
        combinedSequence.addBlock(cBlock)
    else
        for roIdx=1:length(readoutSequence.blockEvents)
            combinedSequence.addBlock(readoutSequence.getBlock(roIdx));
        end
    end
end

%% copy definitions
defs = prepSequence.definitions; % read headers (prep seq)
defKeys = keys(defs);
for ll=1:length(defKeys) % put old headers (readout) in new seq
    combinedSequence.setDefinition(defKeys{ll}, defs(defKeys{ll}));
end

defs = readoutSequence.definitions; % read headers (readout seq)
defKeys = keys(defs);
% Remove 'TotalDuration' from the defKeys cell array
indexToRemove = strcmp(defKeys, 'TotalDuration');
defKeys(indexToRemove) = [];
for ll=1:length(defKeys) 
    combinedSequence.setDefinition(defKeys{ll}, defs(defKeys{ll}));
end

%%
combinedSequence.write(combined_fn);
