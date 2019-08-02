function [StimLength,Omitpre,Omitpost,prefigs,postfigs,red_channel] = SetAnalysisParameters(settingsDir)
%load('C:\Users\bnste\Documents\scripts\Holography_Analysis\Analysis_Parameters.mat');
params = inputdlg({'Enter stimulation length [ms]:','Enter the aquisition frames per second:','Enter # of frames to average pre stimulation:',...
    'Enter # of frames to average post stimulation:' 'Did you record the red channel as well? [Yes-1, No-0]' 'Save as a new default ? [Yes-1, No-0]'}, 'Analysis parameters', 1,...
    {num2str(10),num2str(30),num2str(30),num2str(30),num2str(0),num2str(0)});
StimLength = str2num(params{1});
fps = str2num(params{2});
frametime=1/fps;
Omitpre = 2; % We omit  figures before the stim onset due to stim artifact that might be there
Omitpost = ceil(StimLength*1e-3/frametime)+1;% number of frames to omit after stim. not including stim that is also omitted.
%Omitpost = 3;
prefigs = str2num(params{3});
postfigs = str2num(params{4});
red_channel = str2num(params{5});
if str2num(params{6})==1
    save([settingsDir,'/Analysis_Parameters.mat'],'StimLength','fps','prefigs','postfigs','red_channel', 'Omitpre', 'Omitpost');
end
end