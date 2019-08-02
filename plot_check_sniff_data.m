function [ new_sniff ] = plot_check_sniff_data( m )
%UNTITLED2 Summary of this function goes here
% This function checks the sniff data according to the packet sent time to
% see that no data is missing
samples=m.sniff_samples;time=m.packet_sent_time;sniff=m.sniff;

for idx = 2:size(time)
    checksniff(idx) = time(idx)-time(idx-1)-samples(idx-1);
end

plot(checksniff)
