close all
clear all
clc
figure
d=0.001;
for x=0.05:0.05:0.5
for y=0.05:0.05:0.5
annotation('ellipse',[x,y,d,d],'Facecolor','black');
d=d*1.02;
end
end
