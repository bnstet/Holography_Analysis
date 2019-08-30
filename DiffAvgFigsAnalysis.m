clear variables
close all
clc

[Stack,num_images,fnameStack,fpathStack] = Stack2Figs();%loading the stack of tiff images


[fnamePat, fpathPat]=uigetfile('','Please choose the pattern file','E:\StimTrials');
load([fpathPat fnamePat]);% loading the pattern file

% Seperating the pattern to the different cells and show figs on subplots
[ Xcoordinates , Ycoordinates ] = Cam_spot_builder( pattern, sizesave, xsave ,ysave  );
% for convert=1:size(Xcoordinates,2)
%     [X{convert},Y{convert}] = ds2nfu([(Xcoordinates{convert}-10) Xcoordinates{convert}],[(Ycoordinates{convert}-10) Ycoordinates{convert}]);
% end


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

[ normDiff, Diff ,frameidx,PreStim,PostStim] = DiffAvgFigs( m, Stack ,5 ,15 );

figure
for imgidx = 1:size(frameidx,2)-1
    subplot(ceil((size(frameidx,2)-1)/5),5,imgidx);
    imagesc(Diff(:,:,imgidx));hold on; colormap('gray');
    for convert=1:size(Xcoordinates,2)
        plot(Ycoordinates{convert},Xcoordinates{convert},'or')
    end
    
end

figure;
imagesc(mean(Diff,3));hold on; colormap('gray');
for convert=1:size(Xcoordinates,2)
    plot(Ycoordinates{convert},Xcoordinates{convert},'or')
end

