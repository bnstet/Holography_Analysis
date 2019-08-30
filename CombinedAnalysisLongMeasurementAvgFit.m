clear variables
close all
clc
tic
%% Setting the analysis parameters
[StimLength,Omitpre,Omitpost,prefigs,postfigs] = SetAnalysisParameters();

%% loading the pattern file
[fnamePat, fpathPat]=uigetfile('','Please choose the pattern file','C:\Users\Gilad\Documents\Experiments\Stim trials');
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
[fnameH5, fpathH5]=uigetfile('*.h5','Please choose the hdf5 file of the experiment','C:\Users\Gilad\Documents\Experiments\Stim trials');

% Getting the signals from the HDF5 file, "M" keeps the data in
% Trials structure and "m" concatenates all trials to one colom.
[ M, m ] = HDF5reader( fpathH5,fnameH5,'frame_triggers','sniff','lick1','lick2');
[ M.packet_sent_time, M.sniff_samples, m.packet_sent_time,...
    m.sniff_samples] = HDF5Eventsreader( fpathH5,fnameH5);

%% Checking how many images we have in all the relevant files
[name,path]=uigetfile('*.tif','Please choose the first Tiff Stack files of this session','C:\Users\Gilad\Documents\Experiments\Stim trials');
fname=[path,name];
filenum=1;
while exist(fname,'file')==2
    info = imfinfo(fname);
    num_images(filenum) = numel(info);
    fname(end-4)=num2str(str2num(fname(end-4))+1);
    filenum=filenum+1;
end
totnumimages=sum(num_images);
totnumfiles=filenum-1;
width=info(1).Width; height=info(1).Height;
%% Checking to see that we are not missing data and shift the timing of sniff and figures
[ m.sniff ] = Check_sniff_data( m );
m.time=1:size(m.sniff,1);% This is the local time of the experiment according to sniff recordings
m.ftrig_shift=m.frame_triggers(end-totnumimages+1:end)-single(m.packet_sent_time(1))-single(m.sniff_samples(1));

% Timing the shutter onset
m.shonset_shift=m.shutter_onset-(m.packet_sent_time(1))-(m.sniff_samples(1));% shifting the shutter onset time to sniff time
m.shutter_timing=zeros(1,size(m.sniff,1));
m.shutter_timing(nonzeros(double(m.shonset_shift).*double(m.shonset_shift>0)))=1;

%find closest frame to trigger onset
onsets=find(m.shutter_timing==1);
diffs=min(abs(m.ftrig_shift-find(m.shutter_timing==1)));
for idx1 = 1:size(nonzeros(m.shonset_shift>0),1)
    frameidx(idx1)=find(abs(m.ftrig_shift-onsets(idx1))==diffs(idx1),1); % This is the closest frame to the shutter onset
end
frameidx=frameidx(1:end-1); % remove the last onset

%% open file after file of this session and calculate the relevant parameters

F = zeros(totnumimages,spot_num);
initialFrame=0; % This parameter will hold the frame number between the files
fname=[path,name];
filenum=1;minimum = 0; 
while exist(fname,'file')==2
    [Stack,num_images_temp,fnameStack,fpathStack] = Stack2Figs(name,path);%loading the stack of tiff images
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
    name(end-4)=num2str(str2num(name(end-4))+1);
    fname=[path,name];
    filenum=filenum+1
    initialFrame=frame;
end

%F=F-minimum; % set the lowest value at zero
F=F-min(F(:)); % set the lowest value at zero
    %% Calculating dF/F signal
    F0=prctile(F,5,1);% We define F0 as the lowest 5% of the fluorescence signal over time.
    dF=(F-F0)./F0;

%% Building the plot figure
 Avgimg = mean(Stack,3);
% red = cat(3, ones(size(Stack,1)),zeros(size(Stack,1)), zeros(size(Stack,1)));
% green = cat(3, zeros(size(Stack,1)),ones(size(Stack,1)), zeros(size(Stack,1)));
% 
% transparency = 0.25;
% overlay=zeros(size(pattern,1),size(pattern,2),spot_num);
 h1=figure;
% for idx=1:spot_num
%     img(idx)=subplot(spot_num+1,12,(idx-1)*12+4:(idx-1)*12+5);% Includes saving room for sniff and trig time plot
%     alpha(:,:,idx) = (SpotMat{idx}*transparency);
%     imagesc(Avgimg);colormap('gray');
%     hold on;
%     overlayhandle = imagesc(red);
%     set(overlayhandle, 'AlphaData', alpha(:,:,idx))
%     axis equal;
%     
% end
% img(idx+1)=subplot(spot_num+1,12,(idx)*12+4:(idx)*12+5);imagesc(Avgimg);axis equal;colormap('gray');hold on;
% overlayhandle = imagesc(red);set(overlayhandle, 'AlphaData', sum(alpha,3));
%% Generating and plotting the dF average of all trials per stimulation spot

%[dFvec,frameidx] = AvgFluorescenseTrials( m, dF,0 ); % for global f0
[dFvec,frameidx,Ft0,StimDeltaValue] = AvgFluorescenseTrials( m, F ,1,F0,Omitpre,Omitpost); % for local f0
RealStimVal=0.1;
StimEff=StimDeltaValue>RealStimVal;% Counts all the stims that resulted in df/f> some value
StimAmp=StimDeltaValue.*StimEff;

dFvecmat=cell2mat(dFvec');
for plotidx2=1:size(dFvec,1)
    
    plots2(plotidx2)=subplot(spot_num+1,13,(plotidx2-1)*13+6);
    dFcell = [dFvec{plotidx2,:}];
    dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,:}],2);
    dFSEM(plotidx2,1:size(dFvec{1,1},1))=std(dFcell,0,2)/sqrt(size([dFvec{plotidx2,:}],2));
    ActualStims(plotidx2,1:size(dFvec,2)) = cellfun(@(x) dot(dFavg(plotidx2,:),x)./(norm(dFavg(plotidx2,:))*norm(x)),dFvec(plotidx2,:));
    Tavg=-30+Omitpre:61-Omitpost; % the time for the average df/f plot
    plot(Tavg,dFavg(plotidx2,:),'r','LineWidth',1);
    hold on
    errorshade(Tavg,dFavg(plotidx2,:),dFSEM(plotidx2,:),dFSEM(plotidx2,:),'r','scalar');
    ylim([mean(dFavg(plotidx2,:))-3*std(dFavg(plotidx2,:)) mean(dFavg(plotidx2,:))+3*std(dFavg(plotidx2,:))])
    xlim([Tavg(1)-1 Tavg(end)+1]);
    try
    [fitresult, gof] = FitStim(Tavg, dFavg(plotidx2,:));
    plot(fitresult,'b');legend('off');
    Fits{plotidx2}=fitresult;
    FitCoeffs(plotidx2,:)=coeffvalues(Fits{plotidx2});
    FitAmp(plotidx2)=FitCoeffs(plotidx2,2);
        
    end
    plot([0,0],[get(plots2(plotidx2),'ylim')],'--b')
    ylabel({['Cell # ',num2str(plotidx2)],[num2str(sum(StimEff(plotidx2,:))),'/',...
        num2str(size(StimEff,2)),' stimulated'],['Avg Amp ',num2str(mean(nonzeros(StimAmp(plotidx2,:))))]...
        '\DeltaF/F_0 Avg'},'FontSize',7);
    set(gca,'FontSize',5);
end
FitAmpAvg=mean(FitAmp(FitAmp>0));

%% plotting dF/F , sniff and triggers time
for plotidx=1:spot_num
    MdF=mean(dF(:,plotidx));SdF=std(dF(:,plotidx));
    plots(plotidx)=subplot(spot_num+1,12,(plotidx-1)*12+7:(plotidx-1)*12+12);
%    % d1 = designfilt('bandpassiir','FilterOrder',1, ...
%          'CutoffFrequency1',.1,'CutoffFrequency2',4, ...
%          'SampleRate',30,'DesignMethod','butter');
%    % d1 = designfilt('bandpassfir','FilterOrder',1, ...
%          %'HalfPowerFrequency',0.15,'DesignMethod','butter');
         d1 = designfilt('bandpassfir','FilterOrder',4, ...
         'CutoffFrequency1',0.5,'CutoffFrequency2',4, ...
         'SampleRate',30);
    yfilt = filtfilt(d1,dF(:,plotidx));
    plot(m.ftrig_shift,yfilt,'b');ylim([MdF-2*SdF MdF+2.5*SdF]); ylabel('\DeltaF/F_0');hold on;
%     plot(m.ftrig_shift,dF(:,plotidx),'b');ylim([MdF-2*SdF MdF+2.5*SdF]); ylabel('\DeltaF/F_0');hold on;
    for plotline=1:size(m.time(m.shutter_timing>0),2)
        linetime=m.time(m.shutter_timing>0);
        plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx),'ylim')],'--r')
        if plotline<size(m.time(m.shutter_timing>0),2)
            if StimEff(plotidx,plotline)==1
            rectangle('Position',[linetime(plotline) MdF-2*SdF 500 SdF/5],'EdgeColor','g','FaceColor','g')
            end
        end
    end
end

%plotting the sniff measurements
plots(plotidx+1)=subplot(spot_num+1,12,(plotidx)*12+7:(plotidx)*12+12);
plot(m.time,m.sniff,'b');ylim([min(m.sniff) max(m.sniff)]); ylabel('Sniff');hold on;
for plotline=1:size(m.time(m.shutter_timing>0),2)
    linetime=m.time(m.shutter_timing>0);
    plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx+1),'ylim')],'--r')
end
linkaxes([plots],'x'); % make sure all the x axis limits are the same

sld = uicontrol('Style', 'slider','Units','normalized','Min',0,'Max',1,'SliderStep',[0.01 0.10],'Value',0,'Position', [0.5254    0.0500    0.3796    0.02],'Callback', {@plotxlim,plots,m}); % generate the slider in the plot

%% plotting the difference between pre and post stimulation figures
%figure;
imgtot=subplot(2,12,2:4);
imagesc(Avgimg);hold on; colormap('gray');title('Stack average');


imgmean=subplot(2,12,14:16);
prepost = mean(Diff,3)./max(abs(mean(Diff,3)));
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
toc
