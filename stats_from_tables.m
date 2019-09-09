%% load in experiment stats and info tables

clear all;

tabfile = 'C:\Users\bnste\Downloads\JG1150\JG1150_stat_tables.mat';
load(tabfile);

% movie size
height = 512;
width = 512;

% load red cell indicator list here
redcells = ones(size(maskinds))==0;
redcells(1:52) = true;

% get one-cell-stim recids
onestim_recids = filetab.recid(filetab.onecellstim==1);
%% calculate some useful quantities for later

% get stimmed cells
stimcells = table(unique(pattab.cell), 'VariableNames',{'cell'});


% get table of stimmed cells with the pattern number stimming that cell
% for each recid
patcellmatch = join( removevars(pattab, {'matchpix'}) , stimcells);

% add stimpattern column to stats table (pattern used for that trial)
stattab = join(stattab, removevars(trialtab,{'frame'}));
stattab.Properties.VariableNames{'pattern'} = 'stimpattern';

% add cellpattern column to stats table (pattern corresponding to cell)
stattab = outerjoin( stattab ,patcellmatch, 'RightVariables',{'pattern'});
stattab.Properties.VariableNames{'pattern'} = 'cellpattern';

% add cellstim boolean (whether this cell is stimmed)
stattab.cellstim = (stattab.cellpattern == stattab.stimpattern);

% add redcell boolean (whether this is a red-labeled cell)
stattab.redcell = redcells(stattab.cell);

stattab.onestim = ismember(stattab.recid, onestim_recids);

%% calculate stats for stim cells

cellstim_groups = findgroups(stattab.cell, stattab.cellstim);
cellstim_resp_groups = findgroups(stattab.cell, stattab.cellstim, stattab.resptype);

cellstats = grpstats(stattab, {'cell','recid'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
cellstats_onestim = grpstats(stattab, {'cell', 'onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats = grpstats(stattab, {'cell', 'cellstim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_resptype = grpstats(stattab, {'cell', 'cellstim', 'resptype'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_redcell = grpstats(stattab, { 'cellstim', 'redcell'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_redcell_resptype = grpstats(stattab, { 'cellstim', 'redcell', 'resptype'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});

stimstats_onestim = grpstats(stattab, {'cell', 'cellstim','onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_onestim_resptype = grpstats(stattab, {'cell', 'cellstim', 'resptype','onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});

% get correlation within groups
cellstats.dffstimcorr = splitapply(@(x,y) corr(x,y), double(stattab.respsnr), double(stattab.cellstim), findgroups(stattab.cell, stattab.recid));
cellstats_onestim.dffstimcorr = splitapply(@(x,y) corr(x,y), double(stattab.respsnr), double(stattab.cellstim), findgroups(stattab.cell, stattab.onestim));


tmp = cellstats(~isnan(cellstats.dffstimcorr),:);


%% calculate mean response curves for each cell

mean_resp_stim = {};
std_resp_stim = {};
mean_resp_nostim = {};
std_resp_nostim = {};

ncells=numel(maskinds);

for icell=1:ncells
    tmpstat = stattab(stattab.cell == icell ,{'trial','recid','cell','cellstim'});
    for recind=1:numel(filetab.recid)
        recid = filetab.recid(recind);
        Fc = Fcent_list(recid);
        tmptrialtab = trialtab(trialtab.recid == recid,:);
        inittrial = tmptrialtab.trial(1);
        stimtrialnums = tmpstat.trial(tmpstat.recid == recid &  tmpstat.cellstim);
        nostimtrialnums = tmpstat.trial(tmpstat.recid == recid  & ~tmpstat.cellstim);
        
        tmpstimtrace = Fc(:,stimtrialnums - inittrial + 1, icell);
        tmpnostimtrace = Fc(:,nostimtrialnums - inittrial + 1, icell);
        if numel(tmpstimtrace) > 0
            tmpstimtrace = ( tmpstimtrace - mean(tmpstimtrace(prestiminds,:),1) ) / mean(tmpstimtrace(prestiminds,:),1) ;
        end
        if numel(tmpnostimtrace) > 0
            tmpnostimtrace = ( tmpnostimtrace - mean(tmpnostimtrace(prestiminds,:),1) ) / mean(tmpnostimtrace(prestiminds,:),1) ;
        end     
        if recind==1
            stimtraces = tmpstimtrace;
            nostimtraces = tmpnostimtrace;
        else
            stimtraces = cat(2,stimtraces, tmpstimtrace);
            nostimtraces = cat(2,nostimtraces, tmpnostimtrace);
        end
    end
    mean_resp_stim{icell} = mean(stimtraces,2);
    std_resp_stim{icell} = std(stimtraces,0,2);
    mean_resp_nostim{icell} = mean(nostimtraces,2);
    std_resp_nostim{icell} = std(nostimtraces,0,2);
end


%% plots

% plot stim cell responses
tmp = cat(2,mean_resp_stim{stimcells.cell});
figure; plot( tmp([prestiminds poststiminds],:)); axis([-inf inf -.1 .1]);

tmp = cat(2,mean_resp_nostim{stimcells.cell});
figure; plot( tmp([prestiminds poststiminds],:)); axis([-inf inf -.1 .1]);


%% investigate smoothed trial traces
cellind = 1;
recind = 1;

trialind=1;
window = [1:prestiminds(end) poststiminds(1):size(Fcent_list(1),1)];

%% 

Fc = Fcent_list(recind);

plot(smooth(window, Fc(window,trialind,cellind), 0.3, 'loess')); title(sprintf('movie %d, cell %d, trial %d', recind, cellind, trialind) )
if trialind < size(Fc,2)
    trialind = trialind+1;
else
    trialind = 1;
    cellind = cellind + 1;
end

%% get spatial probability map from multistim experiments


blursigma = 5; % blurring filter

% get relevant vids
msrecs = filetab.recid(filetab.ncellstim == 18 & ~filetab.onecellstim);
multistim_cells = sort(unique(pattab.cell(ismember(pattab.recid,msrecs))));

% get excitation and inhibition probabilities by recid, cell
total_tab = groupsummary(stattab(ismember(stattab.recid, msrecs),{'cell','recid','respsnr'}), {'cell','recid'});
total_counts = total_tab.GroupCount;
excite_tab = groupsummary(stattab(ismember(stattab.recid, msrecs) & stattab.resptype==1,{'cell','recid','respsnr'}), {'cell','recid'}, 'IncludeEmptyGroups',true);
inhibit_tab = groupsummary(stattab(ismember(stattab.recid, msrecs) & stattab.resptype==-1,{'cell','recid','respsnr'}), {'cell','recid'}, 'IncludeEmptyGroups',true);
total_tab.exciteprob = excite_tab.GroupCount ./ total_tab.GroupCount;
total_tab.inhibitprob = inhibit_tab.GroupCount ./ total_tab.GroupCount;


% get mean x,y position for each mask
xavg = zeros(1,numel(maskinds));
yavg = zeros(1,numel(maskinds));
for i=1:numel(maskinds)
    [yinds, xinds] = ind2sub([height width], maskinds{i});
    xavg(i) = mean(xinds(:));
    yavg(i) = mean(yinds(:));
end

% for each mask find the closest stimmed mask and the displacement from
% the stimmed mask to the given mask
closeststim = zeros(1,numel(maskinds));
closestdispx = zeros(1,numel(maskinds));
closestdispy = zeros(1,numel(maskinds));
for i=1:numel(maskinds)
    [~, closeststim(i)] = min( (xavg(multistim_cells) - xavg(i)).^2 + (yavg(multistim_cells) - yavg(i)).^2 );
    closestdispx(i) = round(xavg(i) - xavg(multistim_cells(closeststim(i))));
    closestdispy(i) = round(yavg(i) - yavg(multistim_cells(closeststim(i))));
end




% use displacements to create a weighted mean plot centered on stim sites
newheight = 2*height; newwidth=2*width;
excite_result = zeros(newheight, newwidth);
inhibit_result = zeros(newheight, newwidth);
% go through videos and get stim-centered averages
centery = round(newheight / 2);
centerx = round(newheight / 2);


% create spatial cell weighting image (blurred)
cell_weights = zeros(newheight, newwidth) + 1e-8;
for i=1:numel(maskinds)
    [celly,cellx] =  ind2sub([height width], maskinds{i});
    newx = round(cellx - xavg(i) + centerx + closestdispx(i)); newy = round(celly - yavg(i) + centery +  closestdispy(i));
    newinds = sub2ind([newheight newwidth],newy,newx);
    cell_weights(newinds) = cell_weights(newinds) + 1;
end
cell_weights = imgaussfilt(cell_weights,blursigma);


for irec=1:numel(msrecs)
    recexresult = zeros(newheight, newwidth);
    recinresult = zeros(newheight, newwidth);
    for icell=1:numel(maskinds)
        exweight = total_tab.exciteprob(total_tab.cell == icell & total_tab.recid == msrecs(irec));
        inweight = total_tab.inhibitprob(total_tab.cell == icell & total_tab.recid == msrecs(irec));
        [celly,cellx] =  ind2sub([height width], maskinds{icell});
        newx = round(cellx - xavg(icell) + centerx + closestdispx(icell)); newy = round(celly - yavg(icell) + centery +  closestdispy(icell));
        newinds = sub2ind([newheight newwidth],newy,newx);
        recexresult(newinds) = recexresult(newinds) +  exweight  ;
        recinresult(newinds) = recinresult(newinds) +  inweight   ;
    end
    
    excite_result = excite_result + ( recexresult / numel(msrecs));
    inhibit_result = inhibit_result + ( recinresult / numel(msrecs));
    
end

excite_result = imgaussfilt(excite_result, blursigma) ./ cell_weights;
inhibit_result = imgaussfilt(inhibit_result, blursigma) ./ cell_weights;

diff_result = excite_result - inhibit_result;

climex = .8;
climin = .3;
figure; imagesc(excite_result( (1:height) + centery - round(height/2),(1:width) + centerx - round(width/2))); colorbar; caxis([0 climex]); axis equal;
figure; imagesc(inhibit_result( (1:height) + centery - round(height/2),(1:width) + centerx - round(width/2))); colorbar; caxis([0 climin]); axis equal;
figure; imagesc(diff_result( (1:height) + centery - round(height/2),(1:width) + centerx - round(width/2))); colorbar; caxis([-climin climex]); axis equal
