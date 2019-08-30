function dFvec = stimtrigresponse(m,F,PreStimFrames,PostStimFrames,Omitpre,...
                                           Omitpost,patternTrials,frameidx)
%stimtrigresponse Get the stimulus triggered response for each cell. 
%   Compute dF/F0 for each stimulus and cell.
Ft0 = [];

numCells = size(patternTrials,2);

for cellNum = 1:numCells
    idx2 = 0;
    for idx1 = 1:size(patternTrials{cellNum},1)
        stimNum = patternTrials{cellNum}(idx1);
        if ~ismember(stimNum,m.excludedTrials) % just toss unwanted trials
        idx2 = idx2+1; % pointer to remove zero entries
        stimIdx = find(m.includedTrials==stimNum); % awkward indexing to make up for shift in frameidx
        
        Ft0=mean(F(frameidx(stimIdx)-PreStimFrames:frameidx(stimIdx)-Omitpre,cellNum));% Calculate local f0 for each stimulus
        Fvec{cellNum}{idx2}=[F(frameidx(stimIdx)-PreStimFrames:frameidx(stimIdx)-Omitpre,cellNum);F(frameidx(stimIdx)+Omitpost:frameidx(stimIdx)+PostStimFrames,cellNum)]; % The signal 1 sec before and up to 2 sec after the shutter onset
        dFvec{cellNum}{idx2}=(Fvec{cellNum}{idx2}-Ft0)./Ft0; % Calculate (F-Ft0)/Ft0

        end
    end  
end

%% rotate for ease of handling later
dFvec=dFvec';

