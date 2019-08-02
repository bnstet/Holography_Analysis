function [dFFvec, dFvec , diffStimImages, diffStimLocalImages, localBaselineImages, stimFrames] = stimtrigresponse(m,F,dF,mode,PreStimFrames,PostStimFrames,Omitpre,...
                                           Omitpost,patternTrialsValid,frameidx, Diffidx, ...
                                           PreStim, PostStim, localSize, spots, baselineImg, baselinetype, F0, smoothSize,...
                                           preCalcPeriod, postCalcPeriod)
%stimtrigresponse Get the stimulus triggered response for each cell. 
%   Compute dF/F0 for each cell's activity centered on stimulus times.
% 'baselinetype': 'trial' or 'global'
% mode: 'from_F' or 'from_dF'
halfsize = floor(localSize/2);

Ft0 = [];

numCells = size(patternTrialsValid,2);

for cellNum = 1:numCells
    % local window for averages around stim region
    xc = floor(spots.xcoordsAll(cellNum)); yc = floor(spots.ycoordsAll(cellNum));
    xwindow = [xc - halfsize ,xc - halfsize + localSize - 1 ];
    ywindow = [yc - halfsize , yc - halfsize + localSize - 1];
    stimFrames{cellNum}=[];
    for idx1 = 1:size(patternTrialsValid{cellNum},1)
        stimNum = patternTrialsValid{cellNum}(idx1);
        if  ~isempty(find(stimNum==Diffidx, 1)) % just toss unwanted trials
        stimFrame = frameidx(stimNum);
        stimFrames{cellNum}(end+1)=stimFrame;
        
        preWindow = [max(1,stimFrame-PreStimFrames), max(1,stimFrame - Omitpre)];
        preWindowCalc = [max(1,stimFrame-Omitpre-preCalcPeriod), max(1,stimFrame - Omitpre)];
        postWindow = [min(size(F,1),stimFrame+Omitpost),min(size(F,1),stimFrame+PostStimFrames)];
        postWindowCalc = [min(size(F,1),stimFrame+Omitpost),min(size(F,1),stimFrame+Omitpost+postCalcPeriod)];
        
        
        Ft0=mean(F(preWindowCalc(1):preWindowCalc(2),cellNum));% Calculate local f0 for each stimulus
        
        Fvec{cellNum}{idx1}=[F(preWindow(1):preWindow(2),cellNum); ...
            F(postWindow(1):postWindow(2),cellNum)]; % The signal 1 sec before and up to 2 sec after the shutter onset
        
        dFvec{cellNum}{idx1}=(Fvec{cellNum}{idx1}-Ft0);
        if strcmp(baselinetype,'trial')
            dFFvec{cellNum}{idx1} = dFvec{cellNum}{idx1}./ Ft0;
            diffStimImages{cellNum}{idx1} = imgaussfilt(PostStim(:,:,stimNum) - PreStim(:,:,stimNum), smoothSize)./ imgaussfilt( PreStim(:,:,stimNum), smoothSize);
        elseif strcmp(baselinetype,'global')
            dFFvec{cellNum}{idx1} = dFvec{cellNum}{idx1}./ F0(cellNum);
            diffStimImages{cellNum}{idx1} = imgaussfilt(PostStim(:,:,stimNum) - PreStim(:,:,stimNum), smoothSize)./ baselineImg;
        else
            error('baselinetype must be trial or global');
        end
        
        if strcmp(mode, 'from_dF')
            premean=mean(dF(preWindowCalc(1):preWindowCalc(2),cellNum));
            dFFvec{cellNum}{idx1} = [dF(preWindow(1):preWindow(2),cellNum); ...
                            dF(postWindow(1):postWindow(2),cellNum)] - premean;
        end
        
        % fill in diffStimImages
        
        
        % fill in local difference
        temp = diffStimImages{cellNum}{idx1};
        diffStimLocalImages{cellNum}{idx1} = temp(ywindow(1):ywindow(2),xwindow(1):xwindow(2));
        %diffStimLocalImages{cellNum}{idx1} = PostStim(ywindow(1):ywindow(2),xwindow(1):xwindow(2) ,stimNum) - ...
         %   PreStim(ywindow(1):ywindow(2),xwindow(1):xwindow(2) ,stimNum );
        end
    end  
    
    
    if strcmp(baselinetype,'trial')
        % maybe change this - does it even make sense in this case?
        localBaselineImages{cellNum} = mean(PreStim(ywindow(1):ywindow(2),xwindow(1):xwindow(2),patternTrialsValid{cellNum}),3); 
    else
        localBaselineImages{cellNum} = baselineImg(ywindow(1):ywindow(2),xwindow(1):xwindow(2));
    end
end

%% rotate for ease of handling later
dFvec=dFvec';
dFFvec=dFFvec';
%diffStimImages = diffStimImages';
%diffStimLocalImages = diffStimLocalImages';

