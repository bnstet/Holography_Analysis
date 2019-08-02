function [SLM_pattern,patterns,patternTrials] = HDF5_getPatterns(fpathH5,fnameH5)

H5=h5read([fpathH5,fnameH5],'/Trials');
SLM_pattern = string(H5.SLM_pattern');
SLM_pattern = deblank(SLM_pattern); %deblank gets rid of trailing whitespace
patterns = unique(SLM_pattern);

for idx = 1:length(patterns)
    patternTrials{idx} = find(SLM_pattern==patterns(idx));
end

end