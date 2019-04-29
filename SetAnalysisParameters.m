function [StimLength,Omitpre,Omitpost,prefigs,postfigs,red_channel] = SetAnalysisParameters()
load('/gpfs/home/stetlb01/Holography_Analysis/Analysis_Parameters.mat');
params = inputdlg({'Enter stimulation length [ms]:','Enter the aquisition frames per second:','Enter # of frames to average pre stimulation:',...
    'Enter # of frames to average post stimulation:' 'Did you record the red channel as well? [Yes-1, No-0]' 'Save as a new default ? [Yes-1, No-0]'}, 'Analysis parameters', 1,...
    {num2str(StimLength),num2str(fps),num2str(prefigs),num2str(postfigs),num2str(red_channel),num2str(0)});
StimLength = str2num(params{1});
fps = str2num(params{2});
frametime=1/fps;
Omitpre = 2; % We omit  figures before the stim onset due to stim artifact that might be there
Omitpost = fix(StimLength*1e-3/frametime)+1;% number of frames to omit after stim. not including stim that is also omitted.
prefigs = str2num(params{3});
postfigs = str2num(params{4});
red_channel = str2num(params{5});
if str2num(params{6})==1
    save('C:\Users\Rinberg\Documents\MATLAB\TP stimulation system\Analysis_Parameters.mat','StimLength','fps','prefigs','postfigs','red_channel');
end
end