% CellSeriesAnalysis.m
%   This script plots the responses of neurons targeted individually in the
%   Cell Series experiment. You will need to provide a folder containing
%   holograms of individual neurons, the H5 file specifying stimulation
%   times and target identity for each trial, and the tiffs for a given
%   block of trials. Be sure to specify "red_channel" correctly if you are
%   taking the data right off the rig (vs. after motion correction). 


%tic
%% Setting the analysis parameters
settingsDir = 'C:\Users\bnste\Documents\scripts\Holography_Analysis'; 
loadSettings = 0;
if ~loadSettings
    [StimLength,Omitpre,Omitpost,prefigs,postfigs,red_channel] = ...
        SetAnalysisParameters(settingsDir);
else
    load([settingsDir,'/Analysis_Parameters.mat']);
end

% use 'trial' for faster per-trial baselines, and 'global' for slower
% percentile baseline images
baselinetype = 'trial';


%% loading the pattern file

% FIND A WAY TO ACCUMULATE THE COORDINATES ACROSS HOLOGRAMS FOR
% VISUALIZATION. Make a seperate function. also useful for checking in the
% slm software.

fpathPat=uigetdir('C:\Users\bnste\Documents\scripts\jon_2p_data\','Please choose the folder containing patterns');
fnamesPat = dir(fullfile(fpathPat,'*.mat'));

% Check number of holograms, set up variables
spot_num=length(fnamesPat);
Spotidx{spot_num}=[];
SpotMat=cell(spot_num,1);
SpotMat(:)={zeros(512)};
spots.xcoordsAll = [];
spots.ycoordsAll = [];
spots.sizesAll = [];

width = size(SpotMat{1},1);
height = size(SpotMat{1},2);

% Cycle through each hologram and extract the coordinates. Turn these into
% indices for linearly defining stimulated pixels in images.
for idx1 = 1:spot_num
    
    % Load each pattern in directory
    fnamePat = fnamesPat(idx1).name;
    load(fullfile(fpathPat,fnamePat));% loading the pattern file
    
    spots.xcoordsAll(idx1) = xsave;
    spots.ycoordsAll(idx1) = ysave;
    spots.sizesAll(idx1) = sizesave;
    
    % Get coordinates for each pattern. Coordinates
    [Xcoordinates, Ycoordinates] = ...
        Cam_spot_builder(pattern, sizesave, xsave ,ysave);
    
    spots.XcoordinatesAll{idx1} = Xcoordinates;
    spots.YcoordinatesAll{idx1} = Ycoordinates;
    
    % Test for 1 pixel spots
    [c1,d1]=cellfun(@size ,Xcoordinates);
    [c2,d2]=cellfun(@size ,Ycoordinates);
    
    % if it is 1 pixels spot draws a circle around it for analysis
    if sum(c1.*d1.*c2.*d2)/size(c1,2)==1
        [ Xcoordinates ,Ycoordinates ] = ...
            CircleDrawer( Xcoordinates, Ycoordinates );
    end
    
    
    
    % Create lookup for spots
    Spotidx{idx1}=sub2ind(size(pattern),Xcoordinates{:},...
        Ycoordinates{:});
    SpotMat{idx1}(Spotidx{idx1})=1;
    
    
end
%% loading the H5 file

[fnameH5, fpathH5]=uigetfile('*.h5',...
    'Please choose the hdf5 file of the experiment','C:\Users\bnste\Documents\scripts\jon_2p_data\');

% Extract pharospower? 0 = no, 1 = yes. Don't set to 1 if you don't have
% this field in the data. TODO: check for this instead of toggling.
Pharospower=0;

% Extract the recorded signals from Voyeur file
% *Note* Getting the signals from the HDF5 file, "M" keeps the data in
% Trials structure and "m" concatenates all trials to one column.
[ M, m ] = HDF5reader( fpathH5,fnameH5,Pharospower,...
    'frame_triggers','sniff','lick1','lick2');

% Get the patterns and trials
[SLM_pattern,patterns,patternTrials]= HDF5_getPatterns(fpathH5,fnameH5);

% Get timing of transmission events and packets
[ M.packet_sent_time, M.sniff_samples, m.packet_sent_time,...
    m.sniff_samples] = HDF5Eventsreader( fpathH5,fnameH5);

%% Checking how many images we have in all the relevant files
[name,path]=uigetfile('*.tif',...
    'Please choose the Tiff Stack files of this session','C:\Users\bnste\Documents\scripts\jon_2p_data\');

[tiffInfo]=tiff_info(name,path,red_channel);

%% Checking to see that we are not missing data and shift the timing of sniff and figures
[ m.sniff,m.frame_triggers ] = Check_sniff_data( m );

% This is the local time of the experiment according to sniff recordings
m.time=1:size(m.sniff,1);
m.ftrig_shift=m.frame_triggers(1:end)-single(m.packet_sent_time(1))...
    -single(m.sniff_samples(1));

% Number of trials
m.numTrials = size(m.shutter_onset,1);
m.trial = 1:m.numTrials;

% Timing the shutter onset (shifting the shutter onset time to sniff time)
m.shonset_shift=m.shutter_onset-(m.packet_sent_time(1))...
    -(m.sniff_samples(1));

% Remove stims and trials before the imaging started
m.shutter_timing=zeros(1,size(m.sniff,1));
m.includedTrials = nonzeros(double(m.trial') .* double(m.shonset_shift>0));
m.excludedTrials = setdiff(m.trial,m.includedTrials);
m.shutter_timing(nonzeros(double(m.shonset_shift).*double(m.shonset_shift>0)))=1;

% Find closest frame to trigger onset
onsets=find(m.shutter_timing==1);
diffs=min(abs(m.ftrig_shift-find(m.shutter_timing==1)));
nbad = m.numTrials - size(nonzeros(m.shonset_shift>0),1);

for idx1 = 1:size(diffs,1)
    frameidx(idx1)=find(abs(m.ftrig_shift-onsets(idx1))==diffs(idx1),1); % This is the closest frame to the shutter onset
end

% Remove suspicious trials where the stimulus was not aligned to a frame.
diffThresh = 34; % differences between frame and stim time more than 1 frame.
if any(diffs>diffThresh)
    %frameidx(diffs>diffThresh) = []; % remove stim time
    m.excludedTrials(end+1) = m.includedTrials(diffs>diffThresh); % exclude trial
    m.excludedTrails = sort(m.excludedTrials); % make sure list is ordered
    m.includedTrials(diffs>diffThresh)=[]; % remove from included trials
end

% Trim off the last trial (consider removing this)
m.excludedTrials(end+1) = m.includedTrials(end);
m.includedTrials = m.includedTrials(1:end-1); % remove the last trial
%frameidx=frameidx(1:end-1); % remove the last onset

% remove any trials without full pre-post stim periods
clippedTrialCond = (frameidx < prefigs + 1) |(frameidx > tiffInfo.totnumimages - postfigs);
if any(clippedTrialCond)
    m.excludedTrials= [m.excludedTrials,find(clippedTrialCond)];
    m.excludedTrials = sort(m.excludedTrials);
    m.includedTrials(ismember(m.includedTrials, find(clippedTrialCond))) = [];
end


%create modified patternTrials with only included trial indices
for cellIdx = 1:length(patternTrials)
    cellFrameIdx = patternTrials{cellIdx};
    patternTrialsValid{cellIdx} = cellFrameIdx(ismember(cellFrameIdx, m.includedTrials)) - nbad ;
end


%% Open file after file of this session and calculate the relevant parameters

filterSize = 3; % for DF/F calculations we will apply a square median filter with this width
                   % NOTE: needs to be odd
bkgrndPct = 5; % percentile to use for background determination

initialFrames = cumsum(tiffInfo.num_images) ; initialFrames = [(0);initialFrames];
loopOutput = struct('F',{}, 'normDiff', {},'Diff', {},'PreStim',{}, 'PostStim', {}, ...
        'Diffidx',{}, 'bkgrndImg', {});
for i = 1:length(tiffInfo.filelist)
   loopOutput(i).F = [] ;
end

fname=[path,name];
tic
parfor tiffIdx = 1:length(tiffInfo.filelist)
    disp(['Beginning processing on file ', num2str(tiffIdx),' : ', tiffInfo.filelist(tiffIdx).name]);
    
  [loopOutput(tiffIdx).F, loopOutput(tiffIdx).normDiff, loopOutput(tiffIdx).Diff, ...
     loopOutput(tiffIdx).PreStim, loopOutput(tiffIdx).PostStim, loopOutput(tiffIdx).Diffidx,...
     loopOutput(tiffIdx).bkgrndImg, avgImg(:,:,tiffIdx), medianImg(:,:,tiffIdx),maxImg(:,:,tiffIdx), stdImg(:,:,tiffIdx)] = ...
        HolographyDataLoop(tiffInfo,tiffIdx, frameidx, m,prefigs ,postfigs,Omitpost, red_channel, spot_num, ... 
            Spotidx, SpotMat, baselinetype,filterSize, bkgrndPct);
    
    disp(['Completed processing on file ', num2str(tiffIdx),' : ', tiffInfo.filelist(tiffIdx).name]);

end

% concatenate file results
F = []; normDiff=[]; Diff=[]; PreStim = []; PostStim = []; Diffidx = []; 
for i = 1:length(loopOutput)
    F = cat(1,F,loopOutput(i).F);
    normDiff = cat(3,normDiff, loopOutput(i).normDiff);
    Diff = cat(3,Diff, loopOutput(i).Diff);
    PreStim = cat(3,PreStim, loopOutput(i).PreStim);
    PostStim = cat(3,PostStim, loopOutput(i).PostStim);
    Diffidx = cat(2, Diffidx, loopOutput(i).Diffidx);
end



% combine summary images from individual files
avgImg = mean(avgImg,3);
maxImg = max(maxImg,[],3);
stdImg = sqrt( mean(stdImg.^2,3));
baselineImg = mean(reshape([loopOutput.bkgrndImg],width, height,[]),3);

disp('File loop complete');
toc


clear loopOutput;

%% Calculating dF/F signal
f0window = 5; % average over this many frames to reduce noise in baseline
F0=prctile(movmean(F,f0window,1),5,1);% We define F0 as the lowest 5% of the (window-averaged) fluorescence signal over time.
minF = min(F(:));
dF=(F-F0)./F0;


%% Generating and plotting the dF average of all trials per stimulation spot

% square sidelength for calculating local stim average image (pixels)
localSize = 75;

[dFFvec, dFvec, diffStimImages, diffStimLocalImages, localBaselineImages, stimFrames] = stimtrigresponse( m, F,dF ,prefigs,postfigs,...
    Omitpre,Omitpost,patternTrialsValid,frameidx, Diffidx, PreStim, PostStim, localSize, spots, baselineImg, baselinetype, F0, filterSize);

%totLocalMean = mean(reshape([localMeans{:}], localSize,localSize,[]),3);
%dffLocalMeans = cellfun( @(x,y) x ./ y, localMeans, localBaselineImages, 'UniformOutput', false);
dffLocalMeans = cellfun(@(x)  mean(reshape([x{:}], localSize,localSize,[]),3), diffStimLocalImages, 'UniformOutput',false);
totDffLocalMean = mean(reshape([dffLocalMeans{:}], localSize,localSize,[]),3);


%totMeanDiffStimImage = mean(reshape([meanDiffStimImages{:}], width,height,[]),3);
%dffStimChangeImages = cellfun( @(x) x ./ baselineImg , meanDiffStimImages, 'UniformOutput', false);
dffStimChangeImages = cellfun(@(x)  mean(reshape([x{:}], width, height,[]),3), diffStimImages, 'UniformOutput',false);
totDffStimChangeImage = mean(reshape([dffStimChangeImages{:}], width,height,[]),3);

%% Calculate stim success
% RealStimVal=PreSTD;
% StimEff=StimDeltaValue>RealStimVal;% Counts all the stims that resulted in df/f> some value
% StimAmp=StimDeltaValue.*StimEff;
% StimTries=size(StimAmp,2);

% dFFvecmat=cell2mat(dFFvec');

%% Building the plot figure
h1=figure;
TabNum=floor(spot_num/10)+(rem(spot_num,10)>0)*1;
tabgp = uitabgroup(h1,'Position',[0.4 .05 .55 .95]);
for tabs=1:TabNum
    tab(tabs) = uitab(tabgp,'Title',['Cells ',num2str((tabs-1)*10+1),'-',num2str((tabs-1)*10+10), ' Stim Traces']);
end

for plotidx2=1:size(dFFvec,1)
    axes('parent',tab(ceil(plotidx2/10)));
    %plots2(plotidx2)=subplot(spot_num+1,5,(plotidx2-10*floor((plotidx2-1)/10)-1)*13+6);
    plots2(plotidx2)=subplot(10+1,3,(plotidx2-10*floor((plotidx2-1)/10)-1)*3+1:(plotidx2-10*floor((plotidx2-1)/10)-1)*3+2);
    %plots2(plotidx2)=subplot(10,3,plotidx2);
    dFFcell = [dFFvec{plotidx2}{:}];
    dFFavg{plotidx2}=mean([dFFvec{plotidx2}{:}],2);
    dFFSEM{plotidx2}=std(dFFcell,0,2)/sqrt(size([dFFvec{plotidx2}{:}],2));
    %     ActualStims(plotidx2,1:size(dFFvec,2)) = cellfun(@(x) dot(dFavg(plotidx2,:),x)./(norm(dFavg(plotidx2,:))*norm(x)),dFFvec(plotidx2,:));
    Tavg=-prefigs+Omitpre:postfigs+1-Omitpost; % the time for the average df/f plot
    plot(Tavg,dFFavg{plotidx2},'r','Linewidth',1);
    hold on
    errorshade(Tavg,dFFavg{plotidx2},dFFSEM{plotidx2},dFFSEM{plotidx2},'r','scalar');
    ylim([mean(dFFavg{plotidx2})-3*std(dFFavg{plotidx2}) mean(dFFavg{plotidx2})+3*std(dFFavg{plotidx2})])
    xlim([Tavg(1)-1 Tavg(end)+1]);
    
    %     GoodStim(plotidx2)=sum(StimEff(plotidx2,:));
    %     AmpMean(plotidx2)=mean(nonzeros(StimAmp(plotidx2,:)));
    %     AmpSTD(plotidx2)=std(nonzeros(StimAmp(plotidx2,:)));
    
    plot([0,0],[get(plots2(plotidx2),'ylim')],'--b')
    %     ylabel({['Cell # ',num2str(plotidx2)],[num2str(GoodStim(plotidx2)),'/',...
    %         num2str(size(StimEff,2)),' stimulated'],['Avg Amp ',num2str(AmpMean(plotidx2))]...
    %         '\DeltaF/F_0 Avg'},'FontSize',7);
    ylabel({['Cell # ',num2str(plotidx2)],...
        '\DeltaF/F_0 Avg'},'FontSize',7);
    xlabel('# Frames from stim onset');
    set(gca,'FontSize',5);
end


%% plotting the difference between pre and post stimulation figures
% % WILL NEED TO ADJUST THE WAY THINGS ARE DRAWN. This takes a lot of time to
% % execute for some reason. try to fix.
% %figure;


axes('parent',h1)
imgtot=subplot(2,12,2:4);
imagesc(stdImg);hold on; colormap(imgtot,'gray');title('Stack average');
caxis([0 450]);
axis equal
axis tight

plot_spots(spots)
hold off

stimTime = find(Tavg==0);
avgWindow = [dFFavg{:}];
avgWindow = avgWindow(stimTime+1:stimTime+11,:)';
dFFint = mean(avgWindow,2);
for i = 1:size(SpotMat)
    mask_array(:,:,i) = SpotMat{i};
end
responders_idx = 1:size(SpotMat);

resp_ax = subplot(2,12,14:16);
response_map = Response_Map(dFFint,responders_idx,mask_array);
imagesc(response_map)
colormap(resp_ax,'hot')
% caxis([.995 max(unique(response_map))])
caxis([.991 1.1])
axis equal
axis tight
colorbar
title('\DeltaF/F_0 Avg Pre/Post stim');
hold on
plot_spots(spots)

% % savename=name(1:end-10);
% mkdir([path,'Analysis\']);
% clear h1 Stack tempStack ; % don't save the image it is very heavy
% % save([path,'Analysis\',savename,'.mat'],'GoodStim','AmpMean','AmpSTD','StimTries','StimAmp','StimDeltaValue')
% save([path,'Analysis\',savename,'.mat']);
% toc


%% open additional figure windows

clims = [-.12,.12]; % colorbar limits for df/F0
localtickwidth = 10;
localticklabels = [-fliplr(0:localtickwidth:ceil(localSize/2)), localtickwidth:localtickwidth:ceil(localSize/2) ];
localticks = linspace(1, localSize, numel(localticklabels));

for i = 1:spot_num
   newtab = uitab(tabgp,'Title',['Cell ',num2str(i)]);
   axes('parent',newtab)
   %imgtot=subplot(1,3,1);
   %tmpImg = stdImg;
   %x = spots.xcoordsAll(i); y = spots.ycoordsAll(i);
   %tmpImg(y-3:y+3, x-3:x+3)= max(max(stdImg));
   %imagesc(tmpImg);hold on; colormap(imgtot,'gray');title(['Stim spot ', num2str(i)]); 
   
   
   imgDiffPlot = subplot(2,2,1);
   %imgDiff = meanDiffStimImages{i};
   imgDiff = dffStimChangeImages{i};
   imagesc(imgDiff);hold on;  axis equal; axis tight; colormap(imgDiffPlot); caxis(clims);
   title(['global \DeltaF/F_0 stim average: spot ', num2str(i)]);
   plot_spot(spots,i, 'black')
   colorbar; hold off
   
   localDiffPlot = subplot(2,2,2);
   %localDiff = localMeans{i};
   localDiff = dffLocalMeans{i};
   imagesc(localDiff);hold on; axis equal; axis tight;colormap(localDiffPlot);  caxis(clims);
   title('local \DeltaF/F_0 stim image for this spot');
   set(gca, 'XTick', localticks, 'XTickLabel', localticklabels)
   set(gca, 'YTick', localticks, 'YTickLabel', localticklabels)
   colorbar; hold off
   
   localBasePlt = subplot(2,2,3);
   imagesc(localBaselineImages{i});hold on; axis equal; axis tight;colormap(localBasePlt);  
   title('local baseline (F_0) image for this spot');
   set(gca, 'XTick', localticks, 'XTickLabel', localticklabels)
   set(gca, 'YTick', localticks, 'YTickLabel', localticklabels)
   colorbar; hold off
end



%% plot "all stim" summary images
newtab = uitab(tabgp,'Title','All Stims');
axes('parent',newtab)

imgDiffPlot = subplot(2,2,1);
imgDiff = totDffStimChangeImage;
imagesc(imgDiff);hold on; axis on; axis equal; axis tight; colormap(imgDiffPlot); caxis(clims);
title('global \DeltaF/F_0 stim average: all spots ');
colorbar; hold off

localImgDiffPlot = subplot(2,2,2);
imagesc(totDffLocalMean);hold on;  axis equal; axis tight; colormap(localImgDiffPlot); caxis(clims);
colorbar; 
title('mean local response (all stim locations)');
set(gca, 'XTick', localticks, 'XTickLabel', localticklabels)
set(gca, 'YTick', localticks, 'YTickLabel', localticklabels)
hold off;

localBasePlt = subplot(2,2,3);
imagesc(baselineImg);hold on; axis equal; axis tight;colormap(localDiffPlot);
title('global baseline image');
colorbar; hold off


