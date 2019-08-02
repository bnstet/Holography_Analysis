function [ oimg,num_green_images,name,path ] = loadtiff2Figs(name, path, sFrame, num2read, redchannel)

% Copyright (c) 2012, YoonOh Tak
% All rights reserved.
% Modified by Eftychios Pnevmatikakis, 01/2017.

% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Gwangju Institute of Science and Technology (GIST), Republic of Korea nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

if nargin < 4 || isempty(num2read)
    num2read = Inf;
end

if nargin < 4 || isempty(sFrame)
    sFrame = 1;
end

if nargin < 2
    [name,path]=uigetfile('*.tif','Please choose the Tiff Stack file','E:\StimTrials');
end

if nargin < 3
    redchannel = 0;
end
fname=[path,name];

tStart = tic;
warn_old = warning('off', 'all'); % To ignore unknown TIFF tag.

%% Change directory
path_parent = pwd;
[pathstr, ~, ~] = fileparts(fname);
if ~isempty(pathstr)
    if ~exist(pathstr, 'dir')
        cd(pathstr);
    end
end

%% Open file
file_opening_error_count = 0;
while ~exist('tiff', 'var')
    try
        tiff = Tiff(fname, 'r');
    catch
        file_opening_error_count = file_opening_error_count + 1;
        pause(0.1);
        if file_opening_error_count > 5 % automatically retry to open for 5 times.
            reply = input('Failed to open the file. Do you wish to retry? Y/n: ', 's');
            if isempty(reply) || any(upper(reply) == 'Y')
                file_opening_error_count = 0;
            else
                error(['Failed to open the file ''' fname '''.']);
            end
        end
    end
end

%% Load image information
tfl = 0; % Total frame length
tcl = 1; % Total cell length
while true
    tfl = tfl + 1; % Increase frame count
    
    iinfo(tfl).w       = tiff.getTag('ImageWidth');
    iinfo(tfl).h       = tiff.getTag('ImageLength');
    iinfo(tfl).spp     = tiff.getTag('SamplesPerPixel');
    iinfo(tfl).color   = iinfo(tfl).spp > 2; % Grayscale: 1(real number) or 2(complex number), Color: 3(rgb), 4(rgba), 6(rgb, complex number), or 8(rgba, complex number)
    iinfo(tfl).complex = any(iinfo(tfl).spp == [2 6 8]);

    if tfl > 1
        % If tag information is changed, make a new cell
        if iinfo(tfl-1).w ~= iinfo(tfl).w || ...
            iinfo(tfl-1).h ~= iinfo(tfl).h || ...
            iinfo(tfl-1).spp ~= iinfo(tfl).spp || ...
            iinfo(tfl-1).color ~= iinfo(tfl).color || ...
            iinfo(tfl-1).complex ~= iinfo(tfl).complex
            tcl = tcl + 1; % Increase cell count
            iinfo(tfl).fid = 1; % First frame of this cell
        else
            iinfo(tfl).fid = iinfo(tfl-1).fid + 1;
        end
    else
        iinfo(tfl).fid = 1; % Very first frame of this file
    end
    iinfo(tfl).cid = tcl; % Cell number of this frame
    
    if tiff.lastDirectory(), break; end;
    tiff.nextDirectory();
end

%% Load image data
num_green_images=tfl/(redchannel+1);
T = tfl;
num2read = min(num2read,T-sFrame+1);

if tcl == 1 % simple image (no cell)
    
    for tfl = sFrame:redchannel+1:sFrame+num2read-1
        tiff.setDirectory(tfl);
        temp = tiff.read();
        if tfl == sFrame
            oimg = zeros([size(temp),T/(redchannel+1)],'like',temp);
        end
        
        if iinfo(tfl).complex
            temp = temp(:,:,1:2:end-1,:) + temp(:,:,2:2:end,:)*1i;
        end
       
        if ~iinfo(tfl).color
            oimg(:,:,iinfo(tfl).fid - sFrame + 1) = temp; % Grayscale image
        else
            oimg(:,:,:,iinfo(tfl).fid - sFrame + 1) = temp; % Color image
        end
    end
else % multiple image (multiple cell)
    oimg = cell(T/(redchannel+1), 1);
    for tfl = 1:redchannel+1:tfl
        tiff.setDirectory(tfl);
        temp = tiff.read();
        if iinfo(tfl).complex
            temp = temp(:,:,1:2:end-1,:) + temp(:,:,2:2:end,:)*1i;
        end
        if ~iinfo(tfl).color
            oimg{iinfo(tfl).cid}(:,:,iinfo(tfl).fid) = temp; % Grayscale image
        else
            oimg{iinfo(tfl).cid}(:,:,:,iinfo(tfl).fid) = temp; % Color image
        end
    end
end

%% Close file
tiff.close();
cd(path_parent);
warning(warn_old);

display(sprintf('The file was loaded successfully. Elapsed time : %.3f s.', toc(tStart)));
end
