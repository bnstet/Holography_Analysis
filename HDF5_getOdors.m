function [odorInfo] = HDF5_getOdors(fpathH5,fnameH5)

H5=h5read([fpathH5,fnameH5],'/Trials');
olfa1_odors = string(H5.olfas0x3Aolfa_00x3Aodor');
olfa1_odors = deblank(olfa1_odors); %deblank gets rid of trailing whitespace
odors = unique(olfa1_odors);

for idx = 1:length(odors)
    odorTrials{idx} = find(olfa1_odors==odors(idx));
end

olfa1_conc = double(H5.odorconc);
concentrations = unique(olfa1_conc);

for idx = 1:length(concentrations)
    concTrials{idx} = find(olfa1_conc==concentrations(idx));
end

i = 1;
for idx1 = 1:length(odors)
    for idx2 = 1:length(concentrations)
        odorConcTrials{i} = find((olfa1_odors==odors(idx1)).*(olfa1_conc==concentrations(idx2)));
        odorConcLabels{i} = [odors(idx1) num2str(concentrations(idx2))];
        i = i+1;
    end
end
odorInfo.odors = odors;
odorInfo.concentrations = concentrations;
odorInfo.odorTrials = odorTrials;
odorInfo.concTrials = concTrials;
odorInfo.odorConcTrials = odorConcTrials;
odorInfo.odorConcLabels = odorConcLabels;

end