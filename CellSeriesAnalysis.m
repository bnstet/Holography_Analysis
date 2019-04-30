% CellSeriesAnalysis.m
%   This script plots the responses of neurons targeted individually in the
%   Cell Series experiment. You will need to provide a folder containing
%   holograms of individual neurons, the H5 file specifying stimulation
%   times and target identity for each trial, and the tiffs for a given
%   block of trials. Be sure to specify "redchannel" correctly if you are
%   taking the data right off the rig (vs. after motion correction). 

clear all
%close all
clc
%tic
%% Setting the analysis parameters
[StimLength,Omitpre,Omitpost,prefigs,postfigs,redchannel] = ...
    SetAnalysisParameters();
%Omitpre=2;
%% loading the pattern file

% FIND A WAY TO ACCUMULATE THE COORDINATES ACROSS HOLOGRAMS FOR
% VISUALIZATION. Make a seperate function. also useful for checking in the
% slm software.
%
fpathPat=uigetdir('/gpfs/data/shohamlab/shared_data/jon_2p_data/','Please choose the folder containing patterns');
fnamesPat = dir(fullfile(fpathPat,'*.mat'));

% Check number of holograms, set up variables
spot_num=length(fnamesPat);
Spotidx{spot_num}=[];
SpotMat=cell(spot_num,1);
SpotMat(:)={zeros(512)};
spots.xcoordsAll = [];
spots.ycoordsAll = [];
spots.sizesAll = [];

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
% TODO: WTF, let's just have one function that reads everything in.

[fnameH5, fpathH5]=uigetfile('*.h5',...
    'Please choose the hdf5 file of the experiment','/gpfs/data/shohamlab/shared_data/jon_2p_data/');
% [fnameH5, fpathH5]=uigetfile('*.h5','Please choose the hdf5 file of the experiment');

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
    'Please choose the Tiff Stack files of this session','/gpfs/data/shohamlab/shared_data/jon_2p_data/');

[tiffInfo]=tiff_info(name,path,redchannel);

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
for idx1 = 1:size(nonzeros(m.shonset_shift>0),1)
    frameidx(idx1)=find(abs(m.ftrig_shift-onsets(idx1))==diffs(idx1),1); % This is the closest frame to the shutter onset
end

% Remove suspicious trials where the stimulus was not aligned to a frame.
diffThresh = 34; % differences between frame and stim time more than 1 frame.
if any(diffs>diffThresh)
    frameidx(diffs>diffThresh) = []; % remove stim time
    m.excludedTrials(end+1) = m.includedTrials(diffs>diffThresh); % exclude trial
    m.excludedTrails = sort(m.excludedTrials); % make sure list is ordered
    m.includedTrials(diffs>diffThresh)=[]; % remove from included trials
end

% Trim off the last trial (consider removing this)
m.excludedTrials(end+1) = m.includedTrials(end);
m.includedTrials = m.includedTrials(1:end-1); % remove the last trial
frameidx=frameidx(1:end-1); % remove the last onset

%% Open file after file of this session and calculate the relevant parameters

F = zeros(tiffInfo.totnumimages,spot_num);
initialFrame=0; % This parameter will hold the frame number between the files
fname=[path,name];
filenum=1;minimum = 0;
tic
% for j = 1:numel(block_names)
for idx = 1:length(tiffInfo.filelist)
    
    % Get path and file name (needed for legacy Stack2Figs)
    path = tiffInfo.pathlist{idx}(1:end-length(tiffInfo.filelist(idx).name));
    thisName = tiffInfo.filelist(idx).name;
    
    % Read in the tiff
    [Stack,num_images_temp,fnameStack,fpathStack] = Stack2Figs(thisName,path,redchannel);

    % Calculating the fluorescence signal from each frame for all cells
    for frame=initialFrame+1:initialFrame+num_images_temp
        
        % Grab each frame and pull out the mean values from each spot
        tempStack=Stack(:,:,frame-initialFrame);
        F(frame,:)=cellfun(@(x) mean(tempStack(x)),Spotidx);
        
        % Store the minimum value of the whole measurement movie
        if min(tempStack(:))<minimum
            minimum = min(Stack(:)); 
        end
    end
%     
    
%     % Onsets of stimulation in this stack of the whole movie
%     frameidxtemp=frameidx((frameidx>initialFrame)&...
%         (frameidx<initialFrame+num_images_temp))-initialFrame;
%     
%     % Initialize Diffidx
%     Diffidx=find((frameidx>initialFrame)&(frameidx<initialFrame+num_images_temp));
    
    % Make the average Post - Pre image
%     [ normDiff(1:width,1:height,Diffidx), Diff(1:width,1:height,Diffidx),PreStim(1:width,1:height,Diffidx),...
%         PostStim(1:width,1:height,Diffidx)] = DiffAvgFigsLongMeasurements( m, Stack ,frameidxtemp ,prefigs ,postfigs,Omitpost );
    

%     tempnum=num2str(100000+str2num(fname(end-8:end-4))+1);
%     fname(end-8:end-4)=tempnum(2:end);
%     name(end-8:end-4)=tempnum(2:end);
    filenum=filenum+1
    initialFrame=frame;
end
toc

%% Calculating dF/F signal
F=F-min(F(:)); % set the lowest value at zero
F0=prctile(F,5,1);% We define F0 as the lowest 5% of the fluorescence signal over time.
dF=(F-F0)./F0;

%% Average image
avgImg = mean(Stack,3);
medianImg = median(Stack,3);
maxImg = max(Stack,[],3);
stdImg = std(single(Stack),0,3);

%% Generating and plotting the dF average of all trials per stimulation spot
PreStimFrames = 30; % the number of frames before stim to include in avg
PostStimFrames = 60; % the number of frames after stim to include in avg 
IntPre = 1;% # of frames before stim to sum for the calculation of the delta in dF before and after stim
IntPost = 10;% # of frames after stim to sum for the calculation of the delta in dF before and after stim

dFvec = stimtrigresponse( m, F ,PreStimFrames,PostStimFrames,...
    Omitpre,Omitpost,patternTrials,frameidx);

%% Calculate stim success
% RealStimVal=PreSTD;
% StimEff=StimDeltaValue>RealStimVal;% Counts all the stims that resulted in df/f> some value
% StimAmp=StimDeltaValue.*StimEff;
% StimTries=size(StimAmp,2);

% dFvecmat=cell2mat(dFvec');

%% Building the plot figure
h1=figure;
TabNum=floor(spot_num/10)+(rem(spot_num,10)>0)*1;
tabgp = uitabgroup(h1,'Position',[0.4 .05 .55 .95]);
for tabs=1:TabNum
    tab(tabs) = uitab(tabgp,'Title',['Cells ',num2str((tabs-1)*10+1),'-',num2str((tabs-1)*10+10)]);
end

for plotidx2=1:size(dFvec,1)
    axes('parent',tab(ceil(plotidx2/10)));
    %plots2(plotidx2)=subplot(spot_num+1,5,(plotidx2-10*floor((plotidx2-1)/10)-1)*13+6);
    plots2(plotidx2)=subplot(10+1,5,(plotidx2-10*floor((plotidx2-1)/10)-1)*5+1:(plotidx2-10*floor((plotidx2-1)/10)-1)*5+2);
    %plots2(plotidx2)=subplot(10,3,plotidx2);
    dFcell = [dFvec{plotidx2}{:}];
    dFavg(plotidx2,:)=mean([dFvec{plotidx2}{:}],2);
    dFSEM(plotidx2,:)=std(dFcell,0,2)/sqrt(size([dFvec{plotidx2}{:}],2));
    %     ActualStims(plotidx2,1:size(dFvec,2)) = cellfun(@(x) dot(dFavg(plotidx2,:),x)./(norm(dFavg(plotidx2,:))*norm(x)),dFvec(plotidx2,:));
    Tavg=-PreStimFrames+Omitpre:PostStimFrames+1-Omitpost; % the time for the average df/f plot
    plot(Tavg,dFavg(plotidx2,:),'r','Linewidth',1);
    hold on
    errorshade(Tavg,dFavg(plotidx2,:),dFSEM(plotidx2,:),dFSEM(plotidx2,:),'r','scalar');
    ylim([mean(dFavg(plotidx2,:))-3*std(dFavg(plotidx2,:)) mean(dFavg(plotidx2,:))+3*std(dFavg(plotidx2,:))])
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

%% plotting dF/F , sniff and triggers time
% for plotidx=1:spot_num
%     axes('parent',tab(ceil(plotidx/10)));
%     MdF=mean(dF(:,plotidx));SdF=std(dF(:,plotidx));
%     plots(plotidx)=subplot(10+1,12,(plotidx-10*floor((plotidx-1)/10)-1)*12+3:...
%                                (plotidx-10*floor((plotidx-1)/10)-1)*12+12);
%     %    % d1 = designfilt('bandpassiir','FilterOrder',1, ...
%     %          'CutoffFrequency1',.1,'CutoffFrequency2',4, ...
%     %          'SampleRate',30,'DesignMethod','butter');
%     %    % d1 = designfilt('bandpassfir','FilterOrder',1, ...
%     %          %'HalfPowerFrequency',0.15,'DesignMethod','butter');
%     d1 = designfilt('bandpassfir','FilterOrder',8, ...
%         'CutoffFrequency1',0.5,'CutoffFrequency2',3, ...
%         'SampleRate',30);
%     yfilt = filtfilt(d1,dF(:,plotidx));
%     
%     plot(m.ftrig_shift,yfilt,'b');
%     ylim([MdF-2*SdF MdF+2.5*SdF]);
%         ylabel('\DeltaF/F_0'); 
%         hold on;
%     %     plot(m.ftrig_shift,dF(:,plotidx),'b');ylim([MdF-2*SdF MdF+2.5*SdF]); ylabel('\DeltaF/F_0');hold on;
%     %     for plotline=2:size(m.time(m.shutter_timing>0),2)-1 % This is because we remove the first and last stims
%     % for plotline=2:size(StimEff,2) % This is because we remove the first and last stims
%     for idx1 = 1:size(patternTrials{plotidx},1)
%         stimNum = patternTrials{plotidx}(idx1);
%         if ~ismember(stimNum,m.excludedTrials) % just toss unwanted trials
%             linetime = m.shonset_shift(stimNum); %use real stim num indexing.
%             plot([linetime, linetime],[get(plots(plotidx),'ylim')],'--r')
%         end
%     end
%     
% %         for plotline=1:size(m.time(m.shutter_timing>0),2)
% %             linetime=m.time(m.shutter_timing>0);
% %             plot([linetime(plotline),linetime(plotline)],...
% %                         [get(plots(plotidx),'ylim')],'--r')
% %             %     if plotline<size(m.time(m.shutter_timing>0),2)
% %             %         if StimEff(plotidx,plotline)==1
% %             %             rectangle('Position',[linetime(plotline) MdF-2*SdF 500 SdF/5],'EdgeColor','g','FaceColor','g')
% %             %         end
% %             %     end
% %         end
% %         
% %         for plotline=1:size(m.time(m.shutter_timing>0),2)
% %             linetime=m.time(m.shutter_timing>0);
% %             plot([linetime(plotline),linetime(plotline)],...
% %                         [get(gca,'ylim')],'--r')
% %             %     if plotline<size(m.time(m.shutter_timing>0),2)
% %             %         if StimEff(plotidx,plotline)==1
% %             %             rectangle('Position',[linetime(plotline) MdF-2*SdF 500 SdF/5],'EdgeColor','g','FaceColor','g')
% %             %         end
% %             %     end
% %         end
% end
% 
% % Plotting the sniff measurements
% for sniffidx=1:TabNum
%     axes('parent',tab(sniffidx));
%     plots(plotidx+sniffidx)=subplot(10+1,12,10*12+3:10*12+12);
%     plot(m.time,m.sniff,'b');
%     ylim([min(m.sniff) max(m.sniff)]); 
%     ylabel('Sniff');
%     hold on;
%     for plotline=1:size(m.time(m.shutter_timing>0),2)
%         linetime=m.time(m.shutter_timing>0);
%         plot([linetime(plotline),linetime(plotline)],...
%                     [get(plots(plotidx+1),'ylim')],'--r')
%     end
% end
% 
% % Make sure all the x axis limits are the same
% linkaxes([plots],'x');
% 
% % Generate the slider in the plot
% sld = uicontrol('Style', 'slider','Units','normalized','Min',0,'Max',1,...
%     'SliderStep',[0.01 0.10],'Value',0,'Position',...
%     [0.5254    0.0500    0.3796    0.02],'Callback', {@plotxlim,plots,m}); 
    

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
avgWindow = dFavg(:,stimTime+1:stimTime+11);
dFint = mean(avgWindow,2);
for i = 1:size(SpotMat)
    mask_array(:,:,i) = SpotMat{i};
end
responders_idx = 1:size(SpotMat);

resp_ax = subplot(2,12,14:16);
response_map = Response_Map(dFint,responders_idx,mask_array);
imagesc(response_map)
colormap(resp_ax,'hot')
% caxis([.995 max(unique(response_map))])
caxis([.991 1.1])
axis equal
axis tight
colorbar
hold on
plot_spots(spots)

% % savename=name(1:end-10);
% mkdir([path,'Analysis\']);
% clear h1 Stack tempStack ; % don't save the image it is very heavy
% % save([path,'Analysis\',savename,'.mat'],'GoodStim','AmpMean','AmpSTD','StimTries','StimAmp','StimDeltaValue')
% save([path,'Analysis\',savename,'.mat']);
% toc
