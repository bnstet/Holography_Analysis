clear variables
close all
clc

[Stack,num_images,fnameStack,fpathStack] = Stack2Figs();%loading the stack of tiff images

[fnamePat, fpathPat]=uigetfile('','Please choose the pattern file','E:\StimTrials');
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
overlay=zeros(size(pattern,1),size(pattern,2),spot_num);
h1=figure;
for idx=1:spot_num
    Spotidx{idx}=sub2ind(size(pattern),Xcoordinates{idx},Ycoordinates{idx});
    SpotMat{idx}(Spotidx{idx})=1;
    img(idx)=subplot(spot_num+1,7,(idx-1)*7+1);% Includes saving room for sniff and trig time plot
    overlaytemp=(SpotMat{idx}==0).*Stack(:,:,round(size(Stack,3)/2));
    overlaytemp(overlaytemp==0)=max(max(overlaytemp));
    overlay(:,:,idx)=overlaytemp;
    imagesc(overlay(:,:,idx));axis equal;
end
Patidx=cat(3, SpotMat{:});Patidx=sum(Patidx,3);
totoverlay=(Patidx==0).*(Stack(:,:,round(size(Stack,3)/2)));totoverlay(totoverlay==0)=max(max(totoverlay));
img(idx+1)=subplot(spot_num+1,7,(idx)*7+1);imagesc(totoverlay);axis equal;
%img(idx+2)=subplot(spot_num+2,7,(idx+1)*7+1);imagesc(totoverlay);axis equal;

% Calculating the fluorescence signal from each frame for all cells
F = zeros(num_images,spot_num);
for frame=1:num_images
    tempStack=Stack(:,:,frame);
    F(frame,:)=cellfun(@(x) sum(tempStack(x)),Spotidx);
end
                                                                                                                                                                                                                                                                                                                                                        
% Getting the signals from the HDF5 file, "Measurements" keeps the data in
% Trials strycture and meas_long concatenates all trials to one colom.
try% trying to find the h5 file automatically
Stackdir=dir([fpathStack,fnameStack]);
Stackdate=Stackdir.date;
h5dir=dir([fpathStack,'*.h5*']);
h5dates={h5dir.date};h5names={h5dir.name};
matchfile=find(strcmp(h5dates,Stackdate));
fpathH5=fpathStack;
fnameH5=h5names{matchfile};
catch
[fnameH5, fpathH5]=uigetfile('*.h5','Please choose the hdf5 file of the experiment','E:\StimTrials');
end

[ M, m ] = HDF5reader( fpathH5,fnameH5,'frame_triggers','sniff','lick1','lick2');
[ M.packet_sent_time, M.sniff_samples, m.packet_sent_time,...
    m.sniff_samples] = HDF5Eventsreader( fpathH5,fnameH5);

% Checking to see that we are not missing data and shift the timing of sniff and figures
[ m.sniff ] = Check_sniff_data( m );
m.time=1:size(m.sniff,1);% This is the local time of the experiment according to sniff recordings
m.ftrig_shift=m.frame_triggers(end-num_images+1:end)-single(m.packet_sent_time(1))-single(m.sniff_samples(1));

% Timing the shutter onset 
m.shonset_shift=m.shutter_onset-(m.packet_sent_time(1))-(m.sniff_samples(1));% shifting the shutter onset time to sniff time
m.shutter_timing=zeros(1,size(m.sniff,1));
m.shutter_timing(nonzeros(double(m.shonset_shift).*double(m.shonset_shift>0)))=1;

F0=prctile(F,5,1);% We define F0 as the lowest 5% of the fluorescence signal over time.
dF=(F-F0)./F0;

% plotting dF/F , sniff and triggers time
for plotidx=1:spot_num
    plots(plotidx)=subplot('Position',[1/7*2,img(plotidx).Position(2),4/7,img(plotidx).Position(4)]);
    plot(m.ftrig_shift,dF(:,plotidx),'b');ylim([min(dF(:,plotidx)) max(dF(:,plotidx))]); ylabel('\DeltaF/F_0');hold on;
    for plotline=1:size(m.time(m.shutter_timing>0),2)
        linetime=m.time(m.shutter_timing>0);
        plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx),'ylim')],'--r')
    end
end

%plotting the sniff measurements
plots(plotidx+1)=subplot('Position',[1/7*2,img(plotidx+1).Position(2),4/7,img(plotidx+1).Position(4)]);
plot(m.time,m.sniff,'b');ylim([min(m.sniff) max(m.sniff)]); ylabel('Sniff');hold on;
for plotline=1:size(m.time(m.shutter_timing>0),2)
        linetime=m.time(m.shutter_timing>0);
        plot([linetime(plotline),linetime(plotline)],[get(plots(plotidx+1),'ylim')],'--r')
end
linkaxes([plots],'x'); % make sure all the x axis limits are the same


% Generating and plotting the dF average of all trials per stimulation spot
%[dFvec,frameidx] = AvgFluorescenseTrials( m, dF,0 ); % for global f0
[dFvec,frameidx,Ft0] = AvgFluorescenseTrials( m, F ,1,F0); % for local f0
h2=figure;
for idx2=1:spot_num
   % img2(idx2)=subplot(spot_num+2,7,(idx2-1)*7+1);% Includes saving room for sniff and trig time plot
        img2(idx2)=subplot(spot_num,4,(idx2)*4-1);% 
    imagesc(overlay(:,:,idx2));axis equal;
end
img2(idx2+1)=subplot(1,2,1);imagesc(totoverlay);axis equal;colormap('gray');

dFvecmat=cell2mat(dFvec');
for plotidx2=1:size(dFvec,1)
    %plots2(plotidx2)=subplot('Position',[1/7*2,img2(plotidx2).Position(2),4/7,img2(plotidx2).Position(4)]);
    plots2(plotidx2)=subplot(spot_num,8,(plotidx2)*8);
    for plotidx3 = 1:size(dFvec,2)-1
        plot(-29:58,dFvec{plotidx2,plotidx3},'Color',[0.5 0.5 0.5]);ylabel('\DeltaF/F_0 Avg');
        hold on;
    end
    
    dFavg(plotidx2,1:size(dFvec{1,1},1))=mean([dFvec{plotidx2,:}],2);
    plot(-29:58,dFavg(plotidx2,:),'m','LineWidth',1); ylim([min(dFvecmat(:,plotidx2)) max(dFvecmat(:,plotidx2))]);   ylabel('\DeltaF/F_0 Avg');
    
end
