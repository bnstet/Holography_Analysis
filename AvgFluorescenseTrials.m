function [ dFvec ,frameidx,Ft0 ,StimDeltaValue,PreSTD] = AvgFluorescenseTrials( m, F,Ft0indicator,PreStimFrames,PostStimFrames,Omitpre,Omitpost,IntPre,IntPost)
%UNTITLED3 Summary of this function goes here
%   This function takes the fluorescense signal and average it over all
%   trials before and after the stimulation onset

onsets=find(m.shutter_timing==1);
diffs=min(abs(m.ftrig_shift-find(m.shutter_timing==1)));
for idx1 = 1:size(nonzeros(m.shonset_shift>0),1)
    frameidx(idx1)=find(abs(m.ftrig_shift-onsets(idx1))==diffs(idx1),1); % This is the closest frame to the shutter onset
end

dFvec=cell(size(frameidx,2)-1,size(F,2));
StimDeltaValue=zeros(size(frameidx,2)-1,size(F,2));
DeltafromStimPoint =zeros(size(frameidx,2)-1,size(F,2));

for idx2 = 1:size(F,2)
    for idx3 = 2:size(frameidx,2)-1
        
%         dFvec{idx3,idx2}=[F(frameidx(idx3)-PreStimFrames:frameidx(idx3)-Omitpre,idx2);F(frameidx(idx3)+Omitpost:frameidx(idx3)+PostStimFrames,idx2)]; % The signal 1 sec before and up to 2 sec after the shutter onset
        if Ft0indicator==1
            Ft0(idx3,idx2)=mean(F(max(frameidx(idx3)-PreStimFrames,1):frameidx(idx3)-Omitpre,idx2));% Calculate local f0 for each stimulus
            %Ft0(idx3,idx2)= dF(frameidx(idx3)-Omitpre,idx2);% local f0 as last pixel before stim (optional)
            Fvec{idx3,idx2}=[F(max(frameidx(idx3)-PreStimFrames,1):frameidx(idx3)-Omitpre,idx2);F(frameidx(idx3)+Omitpost:frameidx(idx3)+PostStimFrames,idx2)]; % The signal 1 sec before and up to 2 sec after the shutter onset
            dFvec{idx3,idx2}=(Fvec{idx3,idx2}-Ft0(idx3,idx2))./Ft0(idx3,idx2); % Calculate (F-Ft0)/Ft0
            %dFvec{idx3,idx2}=(Fvec{idx3,idx2}-Ft0(idx3,idx2))./F0(idx2); % Calculate (F-Ft0)/F0.
            PreStimValue=mean(dFvec{idx3,idx2}(max(PreStimFrames-Omitpre-IntPre+1,1):PreStimFrames-Omitpre));% average of PreStimFrames samples before stim
            PostStimValue=mean(dFvec{idx3,idx2}(max(PreStimFrames-Omitpre+1,1):PreStimFrames-Omitpre+1+IntPost));% average of 15 samples after the stim
            StimDeltaValue(idx3,idx2)=PostStimValue-PreStimValue;
            PreSTD(idx3,idx2)=std(dFvec{idx3,idx2}(1:max(PreStimFrames-Omitpre,1)));% std of PreStimFrames samples before stim
        end
    end
    
end

dFvec(1,:)={zeros(size(dFvec{2,1}))}; % fill the first row with zeros 

%% get rid of the first row which is just because we stated from idx3=2 and it's all zeros
dFvec=dFvec(2:end,:);
Ft0=Ft0(2:end,:);
StimDeltaValue=StimDeltaValue(2:end,:);
PreSTD=PreSTD(2:end,:);

%% rotate for ease of handling later
dFvec=dFvec';
StimDeltaValue=StimDeltaValue';
PreSTD=PreSTD';
    

