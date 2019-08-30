clear all
clc

m.sniff_samples=floor(100*rand(100,1));
m.packet_sent_time=cumsum(m.sniff_samples);

 m.packet_sent_time(2:end)=m.packet_sent_time(2:end)+15;
% m.packet_sent_time(30:end)=m.packet_sent_time(30:end)+20;
% m.packet_sent_time(50:end)=m.packet_sent_time(50:end)+10;
 m.packet_sent_time(95:end)=m.packet_sent_time(95:end)+20;
m.packet_sent_time(100)=m.packet_sent_time(100)+15;



m.sniff=floor(100*rand(sum(m.sniff_samples),1));