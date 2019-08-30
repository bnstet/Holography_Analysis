function [cleanStack] = RemoveStimArtifact( Stack ,frameidx,Omitpost)
%UNTITLED3 Summary of this function goes here
%   This function takes the movie stack and remove the stimulation artifact
%   from it
numbers=[];
for stim=1:1:size(frameidx,2) % Added 1 since I removed 1 from frameindex in the main code
    numbers=[numbers,frameidx(stim)-1:frameidx(stim)+Omitpost];

end
   Stack(:,:,numbers)=[]; % Removing 1 figure before stim and Omitpost figres after stim

cleanStack = Stack;
end


