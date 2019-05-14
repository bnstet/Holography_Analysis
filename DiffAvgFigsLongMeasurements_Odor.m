function [ normDiff,Diff ,PreStim,PostStim] = DiffAvgFigsLongMeasurements_Odor(m,Stack,frameidx,preframes,postframes,omitpre,stimframes)
%UNTITLED3 Summary of this function goes here
%   This function takes the fluorescense signal and average it over all
%   trials before and after the stimulation onset

PreStim = zeros(size(Stack,1),size(Stack,2),size(frameidx,2)-2);
PostStim = zeros(size(Stack,1),size(Stack,2),size(frameidx,2)-2);
for stim=1:1:size(frameidx,2) % use the full range since we already exclude trials
    PreStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)-preframes:frameidx(stim)-omitpre),3); % omit 1 frames before stim to reject artefacts
    
    % Test to see if there are enough frames to average after stim. 
    if frameidx(stim)+(stimframes+postframes)< size(Stack,3)
        PostStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)+stimframes:frameidx(stim)+(stimframes+postframes)),3); % omit frames after stim to reject artefacts
    else
        PostStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)+stimframes:end),3); % omit frames after stim to reject artefacts 
    end
end

Diff=PostStim-PreStim;
normDiff=Diff./PreStim;
normDiff(isinf(normDiff))=mean(mean(mean(normDiff(~isinf(normDiff)))));

end
% dFvec=dFvec';


