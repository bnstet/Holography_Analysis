function [ Output, longvar ] = HDF5reader1Photon( fpathH5,fnameH5,varargin )
%UNTITLED Summary of this function goes here
%   This function reads an HDF5 file and get the asked field.
% Written by Gilad Lerman
H5=h5read([fpathH5,fnameH5],'/Trials');
Output.shutter_onset=H5.laserontime; longvar.shutter_onset=H5.laserontime;% for now this is the way to read shutter onset

for fieldnum=1:size(varargin,2)
    Outputfield_raw{1,H5.trialNumber(end)}=[];
    Outputfield={};
    for trialidx=1:H5.trialNumber(end)
        trialtxt=num2str(10000+trialidx);
        Outputfield_raw{trialidx}=h5read([fpathH5,fnameH5],['/Trial',trialtxt(2:end),'/',varargin{fieldnum}]);
        Outputfield{trialidx}={cat(1,Outputfield_raw{trialidx}{:})};
        Output.(varargin{fieldnum}){trialidx}=cat(1,Outputfield{trialidx}{:});
    end
    longvar.(varargin{fieldnum})=cat(1,Output.(varargin{fieldnum}){:});
end

end


