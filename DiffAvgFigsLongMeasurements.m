function [ normDiff,Diff ,PreStim,PostStim] = DiffAvgFigsLongMeasurements( m, Stack ,frameidx ,preframes ,postframes,stimframes )
%UNTITLED3 Summary of this function goes here
%   This function takes the fluorescense signal and average it over all
%   trials before and after the stimulation onset



PreStim = zeros(size(Stack,1),size(Stack,2),size(frameidx,2)-2);
PostStim = zeros(size(Stack,1),size(Stack,2),size(frameidx,2)-2);
for stim=2:1:size(frameidx,2) % Changed to 2, in case not enough frames before first stim
    PreStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)-preframes:frameidx(stim)-2),3); % omit 1 frames before stim to reject artefacts
    try % this is to catcj the cases where there are not enough frame after last stim in file
    PostStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)+stimframes:frameidx(stim)+(stimframes+postframes)),3); % omit frames after stim to reject artefacts
    catch
        if frameidx(stim)+(stimframes+postframes)> size(Stack,3)
                PostStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)+stimframes:end),3); % omit frames after stim to reject artefacts
        end
    end

end

Diff=PostStim-PreStim;
normDiff=Diff./PreStim;
normDiff(isinf(normDiff))=mean(mean(mean(normDiff(~isinf(normDiff)))));

end
% dFvec=dFvec';
    

