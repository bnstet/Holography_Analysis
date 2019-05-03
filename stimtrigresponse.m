function [dFvec , diffStimImages, diffStimLocalImages] = stimtrigresponse(m,dF,PreStimFrames,PostStimFrames,Omitpre,...
                                           Omitpost,patternTrials,frameidx, Diffidx, ...
                                           PreStim, PostStim, localSize, spots)
%stimtrigresponse Get the stimulus triggered response for each cell. 
%   Compute dF/F0 for each cell's activity centered on stimulus times.

halfsize = floor(localSize/2);

Ft0 = [];

numCells = size(patternTrials,2);

for cellNum = 1:numCells
    idx2 = 0;
    for idx1 = 1:size(patternTrials{cellNum},1)
        stimNum = patternTrials{cellNum}(idx1);
        if ~ismember(stimNum,m.excludedTrials)  && ~isempty(find(stimNum==Diffidx, 1)) % just toss unwanted trials
        idx2 = idx2+1; % pointer to remove zero entries
        stimIdx = find(m.includedTrials==stimNum); % awkward indexing to make up for shift in frameidx
        
        Ft0=mean(dF(frameidx(stimIdx)-PreStimFrames:...
                frameidx(stimIdx)-Omitpre,cellNum));% Calculate local f0 for each stimulus
        Fvec{cellNum}{idx2}=[dF(frameidx(stimIdx)-PreStimFrames:...
                frameidx(stimIdx)-Omitpre,cellNum);dF(frameidx(stimIdx)+...
                Omitpost:frameidx(stimIdx)+PostStimFrames,cellNum)]; % The signal 1 sec before and up to 2 sec after the shutter onset
        dFvec{cellNum}{idx2}=(Fvec{cellNum}{idx2}-Ft0); % Calculate (F-Ft0)/Ft0 % REMOVED NORMALIZATION
        
        % fill in diffStimImages
        diffStimImages{cellNum}{idx2} = PostStim(:,:,stimNum) - PreStim(:,:,stimNum);
        
        % fill in local difference
        % FIX TRANSPOSE HERE
        xc = floor(spots.ycoordsAll(cellNum)); yc = floor(spots.xcoordsAll(cellNum));
        diffStimLocalImages{cellNum}{idx2} = PostStim(xc - halfsize:xc - halfsize + localSize - 1, ...
            yc - halfsize:yc - halfsize + localSize - 1,stimNum) - ...
            PreStim(xc - halfsize:xc - halfsize + localSize - 1, yc - halfsize:yc - halfsize + localSize - 1,stimNum);
        end
    end  
end

%% rotate for ease of handling later
dFvec=dFvec';
diffStimImages = diffStimImages';
diffStimLocalImages = diffStimLocalImages';

