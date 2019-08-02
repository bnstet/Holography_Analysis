function [ normDiff,Diff ,frameidx,PreStim,PostStim] = DiffAvgFigs( m, Stack ,preframes ,postframes )
%UNTITLED3 Summary of this function goes here
%   This function takes the fluorescense signal and average it over all
%   trials before and after the stimulation onset

onsets=find(m.shutter_timing==1);
diffs=min(abs(m.ftrig_shift-find(m.shutter_timing==1)));
for idx1 = 1:size(nonzeros(m.shonset_shift>0),1)
    frameidx(idx1)=find(abs(m.ftrig_shift-onsets(idx1))==diffs(idx1),1); % This is the closest frame to the shutter onset
end

PreStim = zeros(size(Stack,1),size(Stack,2),size(frameidx,2)-1);
PostStim = zeros(size(Stack,1),size(Stack,2),size(frameidx,2)-1);
for stim=1:1:size(frameidx,2)-1
    PreStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)-preframes:frameidx(stim)-1),3); % omit 1 frames before stim to reject artefacts
    PostStim(:,:,stim) = mean(Stack(:,:,frameidx(stim)+3:frameidx(stim)+postframes),3); % omit 3 frames after stim to reject artefacts
end

Diff=PostStim-PreStim;
normDiff=Diff./PreStim;
normDiff(isinf(normDiff))=mean(mean(mean(normDiff(~isinf(normDiff)))));












% dFvec=cell(size(frameidx,2)-1,size(Stack,2));
% for idx2 = 1:size(Stack,2)
%     for idx3 = 1:size(frameidx,2)-1
%         
%         
%         dFvec{idx3,idx2}=[Stack(frameidx(idx3)-30:frameidx(idx3)-2,idx2);Stack(frameidx(idx3)+2:frameidx(idx3)+60,idx2)]; % The signal 1 sec before and up to 2 sec after the shutter onset
% 
%     end
%     
%     %This part is for the values of last trial but it can be omitted :
%     try
%         dFvec{idx3+1,idx2}=[dF(frameidx(idx3+1)-30:frameidx(idx3+1),idx2)-2;dF(frameidx(idx3+1)+2:frameidx(idx3+1)+60,idx2)]; 
%         %dFvec{idx3+1,idx2}=dF(frameidx(idx3)-30:frameidx(idx3)+60,idx2);
%     catch
%     
%        % dFvec{idx3+1,idx2}=dF(frameidx(idx3)-30:end,idx2);
%         dFvec{idx3+1,idx2}=[dF(frameidx(idx3+1)-30:frameidx(idx3+1)-2,idx2);dF(frameidx(idx3+1)+2:end,idx2)]; 
%     end
end
% dFvec=dFvec';
    

