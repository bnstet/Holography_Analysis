function [dFFvec, dFvec , diffStimImages, diffStimLocalImages, localBaselineImages] = stimtrigresponse(m,F,dF,PreStimFrames,PostStimFrames,Omitpre,...
                                           Omitpost,patternTrialsValid,frameidx, Diffidx, ...
                                           PreStim, PostStim, localSize, spots, baselineImg, baselinetype, F0)
%stimtrigresponse Get the stimulus triggered response for each cell. 
%   Compute dF/F0 for each cell's activity centered on stimulus times.
% 'baselinetype': 'trial' or 'global'
halfsize = floor(localSize/2);

Ft0 = [];

numCells = size(patternTrialsValid,2);

for cellNum = 1:numCells
    for idx1 = 1:size(patternTrialsValid{cellNum},1)
        stimNum = patternTrialsValid{cellNum}(idx1);
        if  ~isempty(find(stimNum==Diffidx, 1)) % just toss unwanted trials
        stimFrame = frameidx(stimNum);
        preWindow = [max(1,stimFrame-PreStimFrames), max(1,stimFrame - Omitpre)];
        postWindow = [min(size(dF,1),stimFrame+Omitpost),min(size(dF,1),stimFrame+PostStimFrames)];
        Ft0=mean(dF(preWindow(1):preWindow(2),cellNum));% Calculate local f0 for each stimulus
        Fvec{cellNum}{idx1}=[dF(preWindow(1):preWindow(2),cellNum); ...
            dF(postWindow(1):postWindow(2),cellNum)]; % The signal 1 sec before and up to 2 sec after the shutter onset
        dFvec{cellNum}{idx1}=(Fvec{cellNum}{idx1}-Ft0);
        if strcmp(baselinetype,'trial')
            dFFvec{cellNum}{idx1} = dFvec{cellNum}{idx1}./ Ft0;
        else
            dFFvec{cellNum}{idx1} = dFvec{cellNum}{idx1}./ F0;
        end
        
        % fill in diffStimImages
        diffStimImages{cellNum}{idx1} = PostStim(:,:,stimNum) - PreStim(:,:,stimNum);
        
        % fill in local difference
        % FIX TRANSPOSE HERE
        xc = floor(spots.xcoordsAll(cellNum)); yc = floor(spots.ycoordsAll(cellNum));
        xwindow = [xc - halfsize ,xc - halfsize + localSize - 1 ];
        ywindow = [yc - halfsize , yc - halfsize + localSize - 1];
        diffStimLocalImages{cellNum}{idx1} = PostStim(ywindow(1):ywindow(2),xwindow(1):xwindow(2) ,stimNum) - ...
            PreStim(ywindow(1):ywindow(2),xwindow(1):xwindow(2) ,stimNum );
        end
    end  
    localBaselineImages{cellNum} = baselineImg(ywindow(1):ywindow(2),xwindow(1):xwindow(2));
end

%% rotate for ease of handling later
dFvec=dFvec';
dFFvec=dFFvec';
%diffStimImages = diffStimImages';
%diffStimLocalImages = diffStimLocalImages';

