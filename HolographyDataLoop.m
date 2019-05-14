function [F, normDiff, Diff, PreStim, PostStim, Diffidx, bkgrndImg, avgImg, medianImg,maxImg, stdImg] = HolographyDataLoop(tiffInfo, tiffIdx, frameidx, m, prefigs, ...
        postfigs,Omitpost, redchannel, spot_num, Spotidx, SpotMat, filterSize, bkgrndPct)

    % processing called on a single loop of the holographic video analysis
    % returns cell fluorescence traces + stim images corresponding to video
    % with index "tiffIdx" in tiffInfo
 
    timeDownsample = 4; % downsample by this factor when calculating baseline
    
    
    width = size(SpotMat{1},1);
    height = size(SpotMat{1},2);
    
    % Get path and file name (needed for legacy Stack2Figs)
    path = tiffInfo.pathlist{tiffIdx}(1:end-length(tiffInfo.filelist(tiffIdx).name));
    thisName = tiffInfo.filelist(tiffIdx).name;

    % get initial frame index
    initialFrames = cumsum(tiffInfo.num_images) ; initialFrames = [(0);initialFrames];
    initialFrame = initialFrames(tiffIdx);
     
    % Read in the tiff
    [Stack,num_images_temp,fnameStack,fpathStack] = Stack2Figs(thisName,path,redchannel);
    
    
    avgImg = mean(Stack,3);
    medianImg = median(Stack,3);
    maxImg = max(Stack,[],3);
    stdImg = std(single(Stack),0,3);
    
    tic
    % Make the background image based of a median-filtered stack
    tInds = randsample(size(Stack,3), floor(size(Stack,3) / timeDownsample));
    bkgrndImg = prctile(imgaussfilt3(Stack(:,:,tInds),filterSize/2, 'FilterSize', [filterSize,filterSize,1]), bkgrndPct, 3 );
    toc
    
    F = zeros(num_images_temp,spot_num);

    % Calculating the fluorescence signal from each frame for all cells
    for frame=initialFrame+1:initialFrame+num_images_temp
        
        % Grab each frame and pull out the mean values from each spot
        tempStack=Stack(:,:,frame-initialFrame);
        F(frame - initialFrame,:)=cellfun(@(x) mean(tempStack(x)),Spotidx);
        
    end
%     
    
     % Onsets of stimulation in this stack of the whole movie
     frameidxtemp=frameidx((frameidx>initialFrame)&...
         (frameidx<initialFrame+num_images_temp))-initialFrame;
     
     % Initialize Diffidx
    Diffidx=find((frameidx>initialFrame)&(frameidx<initialFrame+num_images_temp));
    
    % Make the average Post - Pre image
     [ normDiff(1:width,1:height,1:length(Diffidx)), Diff(1:width,1:height,1:length(Diffidx)),PreStim(1:width,1:height,1:length(Diffidx)),...
         PostStim(1:width,1:height,1:length(Diffidx))] = DiffAvgFigsLongMeasurements( m, Stack ,frameidxtemp ,prefigs ,postfigs,Omitpost);
     
     
     
    
end

