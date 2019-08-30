clear variables
close all
clc
tic
%% Setting the analysis parameters
[StimLength,Omitpre,Omitpost,prefigs,postfigs,redchannel] = SetAnalysisParameters();
%Omitpre=2;
%% loading the pattern file
[fnamePat, fpathPat]=uigetfile('','Please choose the pattern file','E:\StimTrials');
%[fnamePat, fpathPat]=uigetfile('','Please choose the pattern file','E:\StimTrials\161201\85 KHz');

load([fpathPat fnamePat]);% loading the pattern file
% Seperating the pattern to the different cells and show figs on subplots
[ Xcoordinates , Ycoordinates ] = Cam_spot_builder( pattern, sizesave, xsave ,ysave  );
spot_num=size(Xcoordinates,2);
[c1,d1]=cellfun(@size ,Xcoordinates);[c2,d2]=cellfun(@size ,Ycoordinates);
if sum(c1.*d1.*c2.*d2)/size(c1,2)==1  % if it is 1 pixels spot draws a circle around it for analysis
    [ Xcoordinates ,Ycoordinates ] = CircleDrawer( Xcoordinates, Ycoordinates );
end
Spotidx{spot_num}=[];
SpotMat=cell(spot_num,1); SpotMat(:)={zeros(size(pattern))};
for idx=1:spot_num
    Spotidx{idx}=sub2ind(size(pattern),Xcoordinates{idx},Ycoordinates{idx});
    SpotMat{idx}(Spotidx{idx})=1;
end
%% loading the H5 file
[fnameH5, fpathH5]=uigetfile('*.h5','Please choose the hdf5 file of the experiment','C:\VoyeurData');
%[fnameH5, fpathH5]=uigetfile('*.h5','Please choose the hdf5 file of the experiment','E:\StimTrials\161201\85 KHz');


% Getting the signals from the HDF5 file, "M" keeps the data in
% Trials structure and "m" concatenates all trials to one colom.
Pharospower=0; % if this is one then it will read the pharos power field in h5 file
[ M, m, H5] = HDF5reader( fpathH5,fnameH5,Pharospower,'frame_triggers','sniff','lick1','lick2');
[ M.packet_sent_time, M.sniff_samples, m.packet_sent_time,...
    m.sniff_samples] = HDF5Eventsreader( fpathH5,fnameH5);
[aligned_sniff, all_sniff] = SniffReader( fpathH5,fnameH5 );
%% Checking how many images we have in all the relevant files
[name,path]=uigetfile('*.tif','Please choose the first Tiff Stack files of this session','E:\StimTrials');
%[name,path]=uigetfile('*.tif','Please choose the first Tiff Stack files of this session','E:\StimTrials\161201\85 KHz');

fname=[path,name];
filenum=1;
while exist(fname,'file')==2
    block_names{filenum} = fname;
    info = imfinfo(fname);
    num_images(filenum) = numel(info);
    tempnum=num2str(100000+str2num(fname(end-16:end-12))+1);
    fname(end-16:end-12)=tempnum(2:end);
    filenum=filenum+1;
end
totnumimages=sum(num_images)/(redchannel+1);
totnumfiles=filenum-1;
width=info(1).Width; height=info(1).Height;
%% Checking to see that we are not missing data and shift the timing of sniff and figures
% finalinterval = m.frame_triggers(end)-m.frame_triggers(end-1); % get the last interval between frame triggers (33 or 34)
% framediffs = [33 34];
% startidx = find(framediffs == finalinterval);
%Make an array of infered frame times assuming the last trigger time is
%right
% infered_triggers(1) = m.frame_triggers(end);
% for i = 1:totnumimages-1
% %     if mod(i,2)==0
%     infered_triggers(i+1) = round(infered_triggers(i)-33.38);
% %         infered_triggers(i+1) = infered_triggers(i)-framediffs(startidx);
% % %     else
% %         infered_triggers(i+1) = infered_triggers(i)-framediffs(mod(startidx,2)+1); %This is how you pick the other element in a 2 element matrix
% %     end
% end
% infered_triggers = flip(infered_triggers);

% [ m.sniff,m.frame_triggers ] = Check_sniff_data_edit( m );
[ m.sniff,m.frame_triggers ] = Check_sniff_data( m );



%test if there are any out of order frame packets (timestamps suddenly
%lower number) and push everything back to accomodate
% 
erroridx = find(diff(m.frame_triggers)<33,1);
if ~isempty(erroridx)
    disp('something wrong with trigger packets')
    framediff = abs(m.frame_triggers(erroridx+1)-m.frame_triggers(erroridx));
    if m.frame_triggers(erroridx+2)-m.frame_triggers(erroridx+1)> 33
        framesubtract = framediff + 33;
    else
        framesubtract = framediff + 34;
    end
   m.frame_triggers(1:erroridx) = m.frame_triggers(1:erroridx)-framesubtract;
%    erroridx = find(diff(m.frame_triggers)<33,1);
end
%     

m.time=1:size(m.sniff,1);% This is the local time of the experiment according to sniff recordings
% m.time=1:size(aligned_sniff,1);% This is the local time of the experiment according to sniff recordings

% m.ftrig_shift=m.frame_triggers(end-totnumimages+1:end)-single(m.packet_sent_time(1))-single(m.sniff_samples(1));
m.ftrig_shift=m.frame_triggers(end-totnumimages+1:end)-single(m.packet_sent_time(1))-single(m.sniff_samples(1));



% Timing the shutter onset
m.shonset_shift=m.shutter_onset-(m.packet_sent_time(1))-(m.sniff_samples(1));% shifting the shutter onset time to sniff time
m.shutter_timing=zeros(1,size(m.sniff,1));
m.shutter_timing(nonzeros(double(m.shonset_shift).*double(m.shonset_shift>0)))=1;

%find closest frame to trigger onset
onsets=find(m.shutter_timing==1);
for i = 1:length(onsets)
    diffs(i)=min(abs(m.ftrig_shift-onsets(i)));
end
for idx1 = 1:size(nonzeros(m.shonset_shift>0),1)
    frameidx(idx1)=find(abs(m.ftrig_shift-onsets(idx1))==diffs(idx1),1); % This is the closest frame to the shutter onset
end

frameidx_alltrials = frameidx;
trialnumber = H5.trialNumber;
stimtrials = find(H5.trialtype==0);
nostimtrials = find(H5.trialtype==1);


% Clean up the index of stimulation frames to exclude the first and last if
% there isn't enough time before or after

if totnumimages-frameidx(end)< 76; % we want 2 seconds after for plotting
    frameidx=frameidx(1:end-1); % remove the last onset
    trialnumber = trialnumber(1:end-1);
end

if frameidx(1)< 30 % we want 1 second before stim for plotting and math
    frameidx = frameidx(2:end); %get rid of the first stim. 
    trialnumber = trialnumber(2:end);
end

% %Temporary adjustment to remove missaligned frames
% frameidx = frameidx(6:end); %get rid of the first stim.
% trialnumber = trialnumber(6:end);
% %%%%%%%

stimtrialvec = intersect(trialnumber,stimtrials);
nostimtrialvec = intersect(trialnumber,nostimtrials);

for j = 1:length(stimtrialvec)
    stimtrialidx(j) = find(trialnumber-stimtrialvec(j) == 0);
end

for j = 1:length(nostimtrialvec)
    nostimtrialidx(j) = find(trialnumber-nostimtrialvec(j) == 0);
end
    
%% open file after file of this session and calculate the relevant parameters

F = zeros(totnumimages,spot_num);
initialFrame=0; % This parameter will hold the frame number between the files
fname=[path,name];
filenum=1;minimum = 0;

Stack_allframes = zeros(width,height,totnumimages);

% while exist(fname,'file')==2
for j = 1:numel(block_names)
    %     [Stack,num_images_temp,fnameStack,fpathStack] = Stack2Figs(name,path,redchannel);%loading the stack of tiff images
    Stack = bigread2(block_names{j},1);
    tempnumframes = size(Stack,3)-1;
    startingframe = initialFrame+1;
    Stack_allframes(:,:,startingframe:startingframe+tempnumframes)= Stack;
    num_images_temp = size(Stack,3);
    frameidxtemp=frameidx((frameidx>initialFrame)&(frameidx<initialFrame+num_images_temp))-initialFrame; % Onsets of stimulation in this stack of the whole movie
    %CleanStack=RemoveStimArtifact( Stack ,frameidxtemp,Omitpost);
    %%% KALMAN FILTER HERE
    % Calculating the fluorescence signal from each frame for all cells
    for frame=initialFrame+1:initialFrame+num_images_temp
        tempStack=Stack(:,:,frame-initialFrame);
        F(frame,:)=cellfun(@(x) mean(tempStack(x)),Spotidx);
        
        if min(tempStack(:))<minimum
            minimum = min(Stack(:)); % this will store the inimum value of the whole measurement movie
        end
    end
    
    %% Calculating the difference between pre and post stimulation figures
    %     prefigs=10; postfigs=10;
    %     stimframes = 3; % number of frames the stimulation appeared in
    
    Diffidx=find((frameidx>initialFrame)&(frameidx<initialFrame+num_images_temp));
    
    [ normDiff(1:width,1:height,Diffidx), Diff(1:width,1:height,Diffidx),PreStim(1:width,1:height,Diffidx),...
        PostStim(1:width,1:height,Diffidx)] = DiffAvgFigsLongMeasurements( m, Stack ,frameidxtemp ,prefigs ,postfigs,Omitpost );
    
    %%
    %     tempnum=num2str(100000+str2num(fname(end-16:end-12))+1);
    %     fname(end-16:end-12)=tempnum(2:end);
    %     name(end-16:end-12)=tempnum(2:end);
    %     filenum=filenum+1
    
    initialFrame=frame;
    Stack = [];
    % end
end
%% Subtract Minimum F signal
%F=F-minimum; % set the lowest value at zero
for k = 1:size(F,2)
    F(:,k)=F(:,k)-min(F(:,k)); % set the lowest value at zero
end
%% Calculating dF/F signal
F0=prctile(F,5,1);% We define F0 as the lowest 5% of the fluorescence signal over time.
for i = 1:length(F0)
    dF(:,i)=(F(:,i)-F0(i))./F0(i);
end

%% Building the plot figure
Avgimg = mean(Stack_allframes,3);
h1=figure;
%% Generating and plotting the dF average of all trials per stimulation spot
PreStimFrames = 30; % the number of frames before stim to include in avg dF plots
PostStimFrames = 76; % the number of frames after stim to include in avg dF plots

%[dFvec,frameidx] = AvgFluorescenseTrials( m, dF,0 ); % for global f0
[dFvec,frameidx,Ft0,StimDeltaValue,PreSTD] = AvgFluorescenseTrials( m, F , frameidx,1,PreStimFrames,PostStimFrames,Omitpre,Omitpost); % for local f0
RealStimVal=PreSTD;
StimEff=StimDeltaValue>RealStimVal;% Counts all the stims that resulted in df/f> some value
StimAmp=StimDeltaValue.*StimEff;
StimTries=size(StimAmp,2);

% dFvecmat=cell2mat(dFvec');

TabNum=floor(spot_num/10)+(rem(spot_num,10)>0)*1;
tabgp = uitabgroup(h1,'Position',[0.4 .05 .55 .95]);
for tabs=1:TabNum
    tab(tabs) = uitab(tabgp,'Title',['Cells ',num2str((tabs-1)*10+1),'-',num2str((tabs-1)*10+10)]);
end
% for plotidx2=1:size(dFvec,1)
%     axes('parent',tab(ceil(plotidx2/10)));
%     %plots2(plotidx2)=subplot(spot_num+1,13,(plotidx2-10*floor((plotidx2-1)/10)-1)*13+6);
%     plots2(plotidx2)=subplot(10+1,20,(plotidx2-10*floor((plotidx2-1)/10)-1)*20+1:(plotidx2-10*floor((plotidx2-1)/10)-1)*20+2);
%     %plots2(plotidx2)=subplot(10,3,plotidx2);
%     dFcell = [dFvec{plotidx2,:}];
% %     dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,:}],2);
%     dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,:}],2);
%     dFSEM(plotidx2,1:size(dFvec{1,1},1))=std(dFcell,0,2)/sqrt(size([dFvec{plotidx2,:}],2));
%     ActualStims(plotidx2,1:size(dFvec,2)) = cellfun(@(x) dot(dFavg(plotidx2,:),x)./(norm(dFavg(plotidx2,:))*norm(x)),dFvec(plotidx2,:));
%     Tavg=-PreStimFrames+Omitpre:PostStimFrames+1-Omitpost; % the time for the average df/f plot
%     plot(Tavg,dFavg(plotidx2,:),'r','LineWidth',1);
%     hold on
%     errorshade(Tavg,dFavg(plotidx2,:),dFSEM(plotidx2,:),dFSEM(plotidx2,:),'r','scalar');
%     ylim([mean(dFavg(plotidx2,:))-3*std(dFavg(plotidx2,:)) mean(dFavg(plotidx2,:))+3*std(dFavg(plotidx2,:))])
%     xlim([Tavg(1)-1 Tavg(end)+1]);
%     
%     GoodStim(plotidx2)=sum(StimEff(plotidx2,:));
%     AmpMean(plotidx2)=mean(nonzeros(StimAmp(plotidx2,:)));
%     AmpSTD(plotidx2)=std(nonzeros(StimAmp(plotidx2,:)));
%     
%     plot([0,0],[get(plots2(plotidx2),'ylim')],'--b')
%     %     ylabel({['Cell # ',num2str(plotidx2)],[num2str(GoodStim(plotidx2)),'/',...
%     %         num2str(size(StimEff,2)),' stimulated'],['Avg Amp ',num2str(AmpMean(plotidx2))]...
%     %         '\DeltaF/F_0 Avg'},'FontSize',7);
%     ylabel({['Cell # ',num2str(plotidx2)],...
%         '\DeltaF/F_0 Avg'},'FontSize',7);
%     xlabel('# Frames from stim onset');
%     set(gca,'FontSize',5);
% end   
% d1 = designfilt('bandpassfir','FilterOrder',4, ...
%         'CutoffFrequency1',0.5,'CutoffFrequency2',4, ...
%         'SampleRate',30);

for plotidx2=1:size(dFvec,1)
    axes('parent',tab(ceil(plotidx2/10)));
    %plots2(plotidx2)=subplot(spot_num+1,13,(plotidx2-10*floor((plotidx2-1)/10)-1)*13+6);
    plots2(plotidx2)=subplot(10+1,20,(plotidx2-10*floor((plotidx2-1)/10)-1)*20+1:(plotidx2-10*floor((plotidx2-1)/10)-1)*20+2);
    %plots2(plotidx2)=subplot(10,3,plotidx2);
    dFcell = [dFvec{plotidx2,stimtrialidx}];
    dFcell_nostim = [dFvec{plotidx2,nostimtrialidx}];
%     dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,:}],2);
    dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,stimtrialidx}],2);
    dFavg_nostim(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,nostimtrialidx}],2);
    
    dFSEM(plotidx2,1:size(dFvec{1,1},1))=std(dFcell,0,2)/sqrt(size([dFvec{plotidx2,stimtrialidx}],2));
    dFSEM_nostim(plotidx2,1:size(dFvec{1,1},1))=std(dFcell_nostim,0,2)/sqrt(size([dFvec{plotidx2,nostimtrialidx}],2));
%     ActualStims(plotidx2,1:size(dFvec,2)) = cellfun(@(x) dot(dFavg(plotidx2,:),x)./(norm(dFavg(plotidx2,:))*norm(x)),dFvec(plotidx2,:));
    Tavg=-PreStimFrames+Omitpre:PostStimFrames+1-Omitpost; % the time for the average df/f plot
    Tavg = Tavg/30;
    
%     for plotidx2 = [1 2 3 4 5 19 26]
%     figure
    plot(Tavg,dFavg(plotidx2,:),'r','LineWidth',1);
    hold on
    plot(Tavg,dFavg_nostim(plotidx2,:),'b','LineWidth',1);
    
    errorshade(Tavg,dFavg(plotidx2,:),dFSEM(plotidx2,:),dFSEM(plotidx2,:),'r','scalar');
    errorshade(Tavg,dFavg_nostim(plotidx2,:),dFSEM_nostim(plotidx2,:),dFSEM_nostim(plotidx2,:),'b','scalar');
%     ylim([mean(dFavg(plotidx2,:))-3*std(dFavg(plotidx2,:)) mean(dFavg(plotidx2,:))+3*std(dFavg(plotidx2,:))])
    ylim([mean(dFavg_nostim(plotidx2,:))-3*std(dFavg_nostim(plotidx2,:)) mean(dFavg_nostim(plotidx2,:))+3*std(dFavg_nostim(plotidx2,:))])

    xlim([Tavg(1)-1 Tavg(end)+1]);
    plot([0,0],[get(plots2(plotidx2),'ylim')],'--r')

%     ylim([-0.15 0.35]) %cell 1-4
%         ylim([-0.15 0.8]) %cell 5
%                 ylim([-0.15 0.4]) %cell 19
%                                 ylim([-0.15 0.6]) %cell 26

      xlim([Tavg(1) Tavg(end)]);

%     plot([0,0],[get(plots2(plotidx2),'ylim')],'--b')
%     plot([0,0],[-0.15 0.6],'--r')
%     set(gca, 'box','off')
%     end
    
    GoodStim(plotidx2)=sum(StimEff(plotidx2,:));
    AmpMean(plotidx2)=mean(nonzeros(StimAmp(plotidx2,:)));
    AmpSTD(plotidx2)=std(nonzeros(StimAmp(plotidx2,:)));
    
    
    %     ylabel({['Cell # ',num2str(plotidx2)],[num2str(GoodStim(plotidx2)),'/',...
    %         num2str(size(StimEff,2)),' stimulated'],['Avg Amp ',num2str(AmpMean(plotidx2))]...
    %         '\DeltaF/F_0 Avg'},'FontSize',7);
    ylabel({['Cell # ',num2str(plotidx2)],...
        '\DeltaF/F_0 Avg'},'FontSize',7);
    xlabel('# Frames from stim onset');
    set(gca,'FontSize',5);
end

%% plotting dF/F , sniff and triggers time
for plotidx=1:spot_num
    axes('parent',tab(ceil(plotidx/10)));
    MdF=mean(dF(:,plotidx));SdF=std(dF(:,plotidx));
    plots(plotidx)=subplot(10+1,12,(plotidx-10*floor((plotidx-1)/10)-1)*12+3:(plotidx-10*floor((plotidx-1)/10)-1)*12+12);
    %    % d1 = designfilt('bandpassiir','FilterOrder',1, ...
    %          'CutoffFrequency1',.1,'CutoffFrequency2',4, ...
    %          'SampleRate',30,'DesignMethod','butter');
    %    % d1 = designfilt('bandpassfir','FilterOrder',1, ...
    %          %'HalfPowerFrequency',0.15,'DesignMethod','butter');
    d1 = designfilt('bandpassfir','FilterOrder',4, ...
        'CutoffFrequency1',0.5,'CutoffFrequency2',4, ...
        'SampleRate',30);
    yfilt = filtfilt(d1,dF(:,plotidx));
    plot(m.ftrig_shift,yfilt,'k');ylim([MdF-2*SdF MdF+2.5*SdF]); ylabel('\DeltaF/F_0');hold on;
 

    %     plot(m.ftrig_shift,dF(:,plotidx),'b');ylim([MdF-2*SdF MdF+2.5*SdF]); ylabel('\DeltaF/F_0');hold on;
    %     for plotline=2:size(m.time(m.shutter_timing>0),2)-1 % This is because we remove the first and last stims
%     for plotline=2:size(StimEff,2) % This is because we remove the first and last stims
%         linetime=m.time(m.shutter_timing>0);
%         plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx),'ylim')],'--r')
%         if plotline<size(m.time(m.shutter_timing>0),2)
%             if StimEff(plotidx,plotline)==1
%                 rectangle('Position',[linetime(plotline) MdF-2*SdF 500 SdF/5],'EdgeColor','g','FaceColor','g')
%             end
%         end
%     end
    
    for plotline=1:size(m.time(m.shutter_timing>0),2)
        linetime=m.time(m.shutter_timing>0);
        if H5.trialtype(plotline)==0
            plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx),'ylim')],'--r')
        else
            plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx),'ylim')],'--b')
        end
    end
    
end

%plotting the sniff measurements
for sniffidx=1:TabNum;
    axes('parent',tab(sniffidx));
    plots(plotidx+sniffidx)=subplot(10+1,12,10*12+3:10*12+12);
    plot(m.time,m.sniff,'b');ylim([min(m.sniff) max(m.sniff)]); ylabel('Sniff');hold on;
    for plotline=1:size(m.time(m.shutter_timing>0),2)
        linetime=m.time(m.shutter_timing>0);
        plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx+1),'ylim')],'--r')
    end
end
linkaxes([plots],'x'); % make sure all the x axis limits are the same

sld = uicontrol('Style', 'slider','Units','normalized','Min',0,'Max',1,'SliderStep',[0.01 0.10],'Value',0,'Position', [0.5254    0.0500    0.3796    0.02],'Callback', {@plotxlim,plots,m}); % generate the slider in the plot

%% plotting the difference between pre and post stimulation figures
%figure;
axes('parent',h1)
imgtot=subplot(2,12,2:4);
imagesc(Avgimg);hold on; colormap('gray');title('Stack average');


imgmean=subplot(2,12,14:16);
Diffmean = mean(Diff,3);
maxDiff = max(max(abs(mean(Diff,3))));
prepost = Diffmean./maxDiff;
imagesc(prepost);hold on; colormap('gray');title(['Average ',num2str(postfigs),' figs post stim - ',num2str(prefigs),' figs pre stim ']);

for convert=1:size(Xcoordinates,2)
    %plot(imgmean,xsave(convert),ysave(convert),'or');hold on
    subplot(2,12,14:16);
    rectangle('Position',[xsave(convert)-sizesave(convert)/2,ysave(convert)-sizesave(convert)/2,sizesave(convert),sizesave(convert)],'Curvature',[1,1],'EdgeColor','y');
    text(xsave(convert)+sizesave(convert)/2,ysave(convert)+sizesave(convert)/2,num2str(convert),'FontSize',8,'Color',[1 1 0]);
    %plot(imgtot,xsave(convert),ysave(convert),'or');hold on
    subplot(2,12,2:4);
    rectangle('Position',[xsave(convert)-sizesave(convert)/2,ysave(convert)-sizesave(convert)/2,sizesave(convert),sizesave(convert)],'Curvature',[1,1],'EdgeColor','y');
    text(xsave(convert)+sizesave(convert)/2,ysave(convert)+sizesave(convert)/2,num2str(convert),'FontSize',8,'Color',[1 1 0]);
    
end

% %% Difference Between Stim and No Stim Sniffs
% 
% h2 = figure;
% % trialtimes =  m.time(m.shutter_timing>0);
% stimtimes = H5.laserontime(H5.trialtype == 0);
% nostimtimes = H5.laserontime(H5.trialtype == 1);
% 
% % trialtimes =  aligned_sniff(m.shutter_timing>0);
% % stimtimes = trialtimes(stimtrialvec); nostimtimes = trialtimes(nostimtrialvec);
% % Tavg=-PreStimFrames+Omitpre:PostStimFrames+1-Omitpost;
% alignedtime = linspace(-2000,2000,4001);
% % 
% % for idx = 1:length(stimtimes)
% %     plot(alignedtime, m.sniff(stimtimes(idx)-2000:stimtimes(idx)+2000),'r')
% %     stim_sniff(idx,:) = m.sniff(stimtimes(idx)-2000:stimtimes(idx)+2000);
% %     hold on
% % end
% 
% for idx = 1:length(stimtimes)
% %     plot(alignedtime, all_sniff(stimtimes(idx)-2000:stimtimes(idx)+2000),'r')
%     stim_sniff(idx,:) = all_sniff(stimtimes(idx)-2000:stimtimes(idx)+2000);
%     hold on
% end
% 
% % for idx = 1:length(H5.laserontime)
% %     plot(alignedtime, all_sniff(H5.laserontime(idx)-2000:H5.laserontime(idx)+2000),'r')
% %     stim_sniff(idx,:) = all_sniff(H5.laserontime(idx)-2000:H5.laserontime(idx)+2000);
% %     hold on
% % end
% 
% plot([0 0],[-1000 1000],'--r')
% hold on
% plot(alignedtime,mean(stim_sniff,1),'r','LineWidth',2)
% 
% stim_sniff_SEM = std(single(stim_sniff),0,1)/sqrt(size(stimtimes,1));
% errorshade(alignedtime,mean(stim_sniff,1),stim_sniff_SEM,stim_sniff_SEM,'r','scalar');
% 
% for idx = 1:length(nostimtimes)
% %     plot(alignedtime, m.sniff(nostimtimes(idx)-2000:nostimtimes(idx)+2000),'b')
%     nostim_sniff(idx,:) = all_sniff(nostimtimes(idx)-2000:nostimtimes(idx)+2000);
% end
% 
% % for idx = 1:length(H5.laserontime)
% %     plot(alignedtime, all_sniff(H5.laserontime(idx)-2000:H5.laserontime(idx)+2000),'r')
% %     stim_sniff(idx,:) = all_sniff(H5.laserontime(idx)-2000:H5.laserontime(idx)+2000);
% %     hold on
% % end
% 
% plot(alignedtime,mean(nostim_sniff,1),'b','LineWidth',2)
% nostim_sniff_SEM = std(single(nostim_sniff),0,1)/sqrt(size(nostimtimes,1));
% errorshade(alignedtime,mean(nostim_sniff,1),nostim_sniff_SEM,nostim_sniff_SEM,'b','scalar');


% figure
% plot(alignedtime,abs(mean(stim_sniff,1)-mean(nostim_sniff,1)),'k')


%     plot(m.time(stimtimes(idx)-2000:stimtimes(idx)+2000), m.sniff(stimtimes(idx)-2000:stimtimes(idx)+2000),'r')

%     plot(m.time(stimtimes(idx)-2000:stimtimes(idx)+2000, m.sniff(stimtimes(idx)-2000:stimtimes(idx)+2000),'r')
%     Tavg=-PreStimFrames+Omitpre:PostStimFrames+1-Omitpost; % the time for the average df/f plot
%     plot(Tavg,dFavg(plotidx2,:),'r','LineWidth',1);
%     hold on
%     plot(Tavg,dFavg_nostim(plotidx2,:),'b','LineWidth',1);
%     plot(m.time,m.sniff,'b');
%     ylim([min(m.sniff) max(m.sniff)]); ylabel('Sniff');

% %% Make a summary movie for each condition
% 
% %% Stimulation Trials
% stimtrials = find(H5.trialtype==0);
% beforestimframes = 150;   
% afterstimframes = 150;
% totalframes = beforestimframes+afterstimframes+1;
% Omitpost=1;
% 
% % define the frames with stim
% stimframes = frameidx_alltrials(stimtrials);
% 
% if totnumimages-stimframes(end)<afterstimframes; % we want 200 frames after for plotting
%     stimframes=stimframes(1:end-1); % remove the last onset
% end
% 
% if stimframes(1)< beforestimframes % we want 1 second before stim for plotting and math
%     stimframes = stimframes(2:end); %get rid of the first stim. 
% end
% 
% % cut out all the trials and put them in one matrix
% 
% AvgStimStack = zeros(width,height,totalframes,length(stimframes));
% 
% for idx = 1:length(stimframes)
%     AvgStimStack(:,:,1:totalframes,idx) = Stack_allframes(:,:,stimframes(idx)-beforestimframes:stimframes(idx)+afterstimframes);
% end
%     
% SummaryMovie_Cond0 = mean(AvgStimStack,4);
% SummaryMovie_Cond0 = uint32(SummaryMovie_Cond0);
% 
% saveastiff(SummaryMovie_Cond0,[fname(1:end-4) '_Summary.tif']);
% 
% % % Make Df/f for the movie
% % DFFStack = zeros(width,height,totalframes,length(stimframes));
% % F0frames = AvgStimStack(:,:,1:beforestimframes-1,:);
% % StackF0 = mean(F0frames,3);
% % for i = 1:totalframes
% %     for j = 1:size(AvgStimStack,4)
% %         DFFStack(:,:,i,j) = AvgStimStack(:,:,i,j)-StackF0(:,:,:,j)./StackF0(:,:,:,j);
% %     end
% % end
% % 
% % DFFMovie_Cond0 = mean(DFFStack,4);
% % DFFMovie_Cond0 = int32(DFFMovie_Cond0);
% 
% % saveastiff(DFFMovie_Cond0,[fname(1:end-4) '_DFF.tif']);
% 
% %% No Stimulation Trials %%%%%%
% 
% stimframes = frameidx_alltrials(nostimtrials);
% 
% if totnumimages-stimframes(end)<afterstimframes; % we want 200 frames after for plotting
%     stimframes=stimframes(1:end-1); % remove the last onset
% end
% 
% if stimframes(1)< beforestimframes % we want 1 second before stim for plotting and math
%     stimframes = stimframes(2:end); %get rid of the first stim. 
% end
% 
% % cut out all the trials and put them in one matrix
% 
% AvgStimStack = zeros(width,height,totalframes,length(stimframes));
% 
% for idx = 1:length(stimframes)
%     AvgStimStack(:,:,1:totalframes,idx) = Stack_allframes(:,:,stimframes(idx)-beforestimframes:stimframes(idx)+afterstimframes);
% end
%     
% SummaryMovie_Cond1 = mean(AvgStimStack,4);
% SummaryMovie_Cond1 = uint32(SummaryMovie_Cond1);
% 
% saveastiff(SummaryMovie_Cond1,[fname(1:end-4) '_NoStimSummary.tif']);
% 
% % Make Df/f for the movie
% % DFFStack = zeros(width,height,totalframes,length(stimframes));
% % F0frames = AvgStimStack(:,:,1:beforestimframes-1,:);
% % StackF0 = mean(F0frames,3);
% % for i = 1:totalframes
% %     for j = 1:size(AvgStimStack,4)
% %         DFFStack(:,:,i,j) = AvgStimStack(:,:,i,j)-StackF0(:,:,:,j)./StackF0(:,:,:,j);
% %     end
% % end
% % 
% % DFFMovie_Cond1 = mean(DFFStack,4);
% % DFFMovie_Cond1 = int32(DFFMovie_Cond1);
% % 
% % saveastiff(DFFMovie_Cond1,[fname(1:end-4) '_NoStimDFF.tif']);
% 
% 
% 
% %%
% % savename=name(1:end-10);
% % mkdir([path,'Analysis\']);
% % save([path,'Analysis\',savename,'.mat'],'GoodStim','AmpMean','AmpSTD','StimTries','StimAmp','StimDeltaValue')
% toc
