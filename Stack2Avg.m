function [ AvgFig ] = Stack2Avg(name,path )
%UNTITLED Summary of this function goes here
%   This function gets a stack of alternating tiff figures of the red and green
%   channels and averages the two channels and output an averaged RGB figure. 
if nargin~=2
    [name,path]=uigetfile;
end

fname=[path,name];
info = imfinfo(fname);
num_images = numel(info);
A=zeros(info(1).Width,info(1).Height,num_images);
for k = 1:num_images
    A(:,:,k) = imread(fname, k, 'Info', info);
end

Green=A(:,:,1:2:end);
Red=A(:,:,2:2:end);
AvgG=mean(Green,3);
AvgR=mean(Red,3);

figure; imagesc(AvgG);title('Green channle');
figure;imagesc(AvgR);title('Red channle');

Gn = AvgG-min(min(AvgG));
Gn = Gn./max(max(Gn));

Rn = AvgR-min(min(AvgR));
Rn = Rn./max(max(Rn));

AvgFig=cat(3,Rn,Gn,zeros(info(1).Width,info(1).Height));

figure;imagesc(AvgFig);title('Avg Image - Red + Green');

imwrite(AvgFig,[fname(1:end-4),'Avg.tif']);

end

