function [ new_sniff, triggers ] = Check_sniff_data( m )
%UNTITLED2 Summary of this function goes here
% This function checks the sniff data according to the packet sent time to
% see that no data is missing
samples=m.sniff_samples;time=m.packet_sent_time;sniff=m.sniff;

new_sniff=zeros(time(end)-time(1)+samples(1),1);
new_sniff(1:samples(1))=sniff(1:samples(1));
delta=time-samples;
idx2=1;idx3=1;misses=0;
for idx1=2:size(time)
    prev=time(idx1-1);
    if delta(idx1)==prev
        new_sniff(samples(1)+idx2:samples(1)+idx2+samples(idx1)-1)=...
            sniff(samples(1)+idx3:samples(1)+idx3+samples(idx1)-1);
        idx2=idx2+samples(idx1);
        idx3=idx3+samples(idx1);
    else
        misses=misses+1;
        d=delta(idx1)-prev;
        new_sniff(samples(1)+idx2+d:samples(1)+idx2+d+samples(idx1)-1)=...
            sniff(samples(1)+idx3:samples(1)+idx3+samples(idx1)-1);
        idx2=idx2+d+samples(idx1);
        idx3=idx3+samples(idx1);
    end
end

if misses>0
    e1= msgbox(['There were ',num2str(misses),' missing sniff packets'],...
                                                        'Error','error');
end


triggers=m.frame_triggers;

%% correcting for frame triggers that were written on top of previous frame
%  triggers

% we find the negative values of frame triggers diff and the pos values before them and correct
ftrig_diff = diff(m.frame_triggers);
frameidxneg = find(ftrig_diff<0,1);
while ~isempty(frameidxneg)
    frameidxnegpos = find(ftrig_diff(1:frameidxneg)>35,1,'last');
    fillgap=frameidxneg-frameidxnegpos;
    errorvals = triggers(frameidxnegpos:frameidxneg+1);
    triggers=[triggers(1:frameidxnegpos);...
        round(triggers(frameidxnegpos)+33.3:33.3:...
        triggers(frameidxnegpos)+33.4*fillgap)';...
        triggers(frameidxneg+1:end)]; % make sure to add the right number of frames
    ftrig_diff = diff(triggers);
    fprintf(['We replaced ',num2str(errorvals'),' with ',...
        num2str(triggers(frameidxnegpos:frameidxneg+1)'),'\n']);
    frameidxneg = find(ftrig_diff<0,1);
end

%% correcting for missing frame triggers because of missing sent packets

if misses>0 % If there is a problem with sent packets then we correct the frame triggers
    frameidx=find(diff(triggers)>34,1);
    while ~isempty(frameidx)
        fill=round((triggers(frameidx+1)-triggers(frameidx))/33.3);
        errorvals = triggers(frameidx:frameidx+1);
        triggers=[triggers(1:frameidx);...
            round(triggers(frameidx)+33.3:33.3:triggers(frameidx)+33.4*(fill-1))';...
            triggers((frameidx+1):end)]; % make sure to add the right number of frames
        fprintf(['We replaced ',num2str(errorvals'),' with ',...
            num2str(triggers(frameidx:frameidx+fill)'),'\n']);
        frameidx=find(diff(triggers)>35,1);
    end
end























