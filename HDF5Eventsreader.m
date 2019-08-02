function [ packet_sent_time, sniff_samples,longvar_time,longvar_sniff_samples] = HDF5Eventsreader( fpathH5,fnameH5 )
%UNTITLED Summary of this function goes here
%   This function reads an HDF5 file and get the "Events" field of the file.
% Written by Gilad Lerman
H5=h5read([fpathH5,fnameH5],'/Trials');

    Outputfield_raw{1,H5.trialNumber(end)}=[];
    Outputfield={};
for trialidx=1:H5.trialNumber(end)
    trialtxt=num2str(10000+trialidx);
    Outputfield_raw{trialidx}=h5read([fpathH5,fnameH5],['/Trial',trialtxt(2:end),'/','Events']);
    Outputfield{trialidx}={cat(1,Outputfield_raw{trialidx}.packet_sent_time)};
    packet_sent_time{trialidx}=cat(1,Outputfield{trialidx}{:});
    Outputfield{trialidx}={cat(1,Outputfield_raw{trialidx}.sniff_samples)};
    sniff_samples{trialidx}=cat(1,Outputfield{trialidx}{:});
end
longvar_time=cat(1,packet_sent_time{:});
longvar_sniff_samples=cat(1,sniff_samples{:});
end
