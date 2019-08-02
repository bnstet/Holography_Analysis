function [ Figs,num_green_images,name,path ] = Stack2Figs( name, path, redchannel)
%UNTITLED2 Summary of this function goes here
%   This function gets a stack of alternating tiff figures of the red and green
%   channels and output a matrix containing all the figures.

if nargin < 2
    [name,path]=uigetfile('*.tif','Please choose the Tiff Stack file','E:\StimTrials');
end

if nargin < 3
    redchannel = 0;
end

fname=[path,name];
info = imfinfo(fname);
num_images = numel(info);
num_green_images=num_images/(redchannel+1);
Figs=zeros(info(1).Width,info(1).Height,num_green_images,'int16');
n=1;
for k = 1:redchannel+1:num_images
    Figs(:,:,n) = imread(fname, k, 'Info', info);
    n=n+1;
end
end

