%% experiment trace file should be loaded as 'expfile'
% starting trial index recorded as 'startTrialIdx'. this will  be the
% offset for the trial indices


% if ~exist('plotOn')
%    plotOn = false;
% end

frameidx = expfile.frameidx;
F = expfile.F;


if ~exist('startTrialIdx')
    startTrialIdx = 1;
end


% match stim sites in SpotMat, Spotidx (and used in patternTrials) to masks
maskIdx=zeros(size(expfile.Spotidx));
maskMatchMeas = zeros(size(expfile.Spotidx));
nOverlaps = {};
for i=1:numel(maskIdx)
    thisSpot = expfile.SpotMat{i}; % this is the stim spot
    nOverlaps{i} = cellfun(@(x) sum(thisSpot(x)), expfile.maskinds); %expfile.maskinds is a cell array of the  linear index lists for the cell masks
    [maskMatchMeas(i) , maskIdx(i)] = max(nOverlaps{i});
end
maskIdx = uint16(maskIdx);


% need to handle non-unique matches - make sure we get unique maskIdx
[uniqueM i j] = unique(maskIdx,'first');
indexToDupes = find(not(ismember(1:numel(maskIdx),i)));
tryCount=0;
if numel(indexToDupes) > 0
    dupe_vals = unique(maskIdx(indexToDupes));
    disp(sprintf("Warning: overlapping stim sites / masks in mouse %s, date %s", expfile.mouse, expfile.date));
    for z=1:numel(dupe_vals)
       disp(sprintf("Multiple mask assignments to mask %d",dupe_vals(z)));
    end
end
while numel(indexToDupes) > 0 && tryCount < 1000
    dupe_vals = unique(maskIdx(indexToDupes));
    for dval_ind=1:numel(dupe_vals)
        dval = dupe_vals(dval_ind);
        val_inds = find(maskIdx == dval);
        [~, bestIdx] = max(maskMatchMeas(val_inds));
        val_inds(bestIdx) = [];
        avail_inds = find(~ismember(1:numel(nOverlaps{1}),maskIdx));
        for kind=1:numel(val_inds)
            k = val_inds(kind);
            laps = nOverlaps{k};
            [v,vinds] = sort(laps, 'descend');
            vinds = vinds(ismember(vinds,avail_inds)); v = v(ismember(vinds,avail_inds));
            maskIdx(k) = vinds(1);
            maskMatchMeas(k) = v(1);
        end
    end
    [uniqueM i j] = unique(maskIdx,'first');
    indexToDupes = find(not(ismember(1:numel(maskIdxTest),i)));
    tryCount = tryCount + 1;
end
if numel(indexToDupes) > 0
    error("Cannot create 1-1 matching of stim sites and masks!");
end


% create vector recording the stim mask number for each
% trial

if isfield(expfile,'patternTrials')
    oneCellStim = true;
else
    oneCellStim = false;
end

if oneCellStim
    [patternTrialsUnroll,~] = find(cell2mat(cellfun(@(x) ismember(1:size(frameidx,2),x), expfile.patternTrials, 'UniformOutput',  false)'));
%     maskTrials = maskIdx(patternTrialsUnroll);
else
    patternTrialsUnroll = ones(1, numel(frameidx));
%     maskTrials = repmat(maskIdx,[1 numel(frameidx)]);
end
    patternTrialsUnroll = uint16(patternTrialsUnroll);
%     maskTrials = uint16(maskTrials);


%% reshape fluorescence traces into windows around stim period

ntrialsorig = size(frameidx,2);
nframes = size(F,1);

fullwindow = [(-fullWindowPreSize:-1) (0:fullWindowPostSize)]' + frameidx;

% eliminate bad trials (window overlaps start/end of experiment)
movlen = size(F,1);
badtrials = any( (fullwindow <= 0) | (fullwindow > movlen), 1);
frameidx(badtrials) = [];
patternTrialsUnroll(badtrials) = [];
fullwindow = [(-fullWindowPreSize:-1) (0:fullWindowPostSize)]' + frameidx; % redo window creation

ntrialsgood = size(frameidx,2);

%prestimwindow = (-preCalcPeriod:-Omitpre)' + frameidx; prestimlen = size(prestimwindow,1);
%poststimwindow = (Omitpost+1:postCalcPeriod)' + frameidx; poststimlen = size(poststimwindow,1);
%fullwindow = cat(1, prestimwindow,poststimwindow);
%prestimwindow = reshape(prestimwindow,1,[]);
%poststimwindow = reshape(poststimwindow,1,[]);
fullwindow = reshape(fullwindow,1,[]);
fullwindow = min(fullwindow,size(F,1));

Fcentered = reshape(F(fullwindow,:),fullWindowPreSize + fullWindowPostSize + 1,[],size(F,2));

% get indices of pre-post periods in the new stim-centered array
prestiminds = fullWindowPreSize + (-preCalcPeriod:-Omitpre);
poststiminds = fullWindowPreSize + 1 + (Omitpost:postCalcPeriod);



%% get stim-centered stats (dff, response snr, response type,..)
prestimmean = single(mean(reshape(Fcentered(prestiminds,:,:),numel(prestiminds),[],size(F,2)),1));
poststimmean = single(mean(reshape(Fcentered(poststiminds,:,:),numel(poststiminds),[],size(F,2)),1));
prestimstd = single(std(reshape(Fcentered(prestiminds,:,:),numel(prestiminds),[],size(F,2)),1));
poststimstd = single(std(reshape(Fcentered(poststiminds,:,:),numel(poststiminds),[],size(F,2)),1));
tscore = single(  (poststimmean - prestimmean) ./ sqrt( (prestimstd.^2 ./ size(prestiminds,2)) + (poststimstd.^2 ./ size(poststiminds,2)) ));


dff = single(squeeze((poststimmean - prestimmean) ./ prestimmean));
respsnr = single(squeeze( (poststimmean - prestimmean) ./ prestimstd)); % dF / stddev of prestim measurement period
resptype = int8(reshape( discretize(respsnr(:), [-inf -1 1 inf]) - 2, size(respsnr))); % -1 for inhibition, 1 for response, 0 otherwise
dffmean = squeeze(mean(dff,1));


%% create stat and trial info tables

[triallabel, celllabel] = ind2sub(size(dff), 1:numel(dff));

trialstattab = table( -1 + startTrialIdx + uint32(triallabel(:)), uint16(celllabel(:)),  dff(:), respsnr(:), resptype(:), prestimmean(:), poststimmean(:), prestimstd(:), poststimstd(:), tscore(:),'VariableNames', {'trial','cell','dff','respsnr','resptype', 'premean','postmean','prestd','poststd', 'tscore'} ); %, 'VariableTypes',{'uint16','uint16','single','single','int8'});
trialinfotab = table( -1 + startTrialIdx + uint32((1:size(dff,1))'), patternTrialsUnroll(:), uint32(frameidx(:)),'VariableNames', {'trial','pattern', 'frame'});
if oneCellStim
    patterninfotab = table( uint16(1:numel(maskIdx))', uint16(maskIdx)', uint32(maskMatchMeas(:)), 'VariableNames',{'pattern','cell', 'matchpix'});
else
    patterninfotab = table( uint16(ones(numel(maskIdx),1)), uint16(maskIdx(:)) , uint32(maskMatchMeas(:)), 'VariableNames',{'pattern','cell', 'matchpix'});
end
%%


