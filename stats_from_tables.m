%% load in experiment stats and info tables

clear all;

tabfile = 'C:\Users\bnste\Downloads\JG1150\JG1150_stat_tables.mat';
load(tabfile);

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
stimstats = grpstats(stattab, {'cell', 'cellstim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_resptype = grpstats(stattab, {'cell', 'cellstim', 'resptype'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_redcell = grpstats(stattab, { 'cellstim', 'redcell'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_redcell_resptype = grpstats(stattab, { 'cellstim', 'redcell', 'resptype'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});

stimstats_onestim = grpstats(stattab_onestim, {'cell', 'cellstim','onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});
stimstats_onestim_resptype = grpstats(stattab_onestim, {'cell', 'cellstim', 'resptype','onestim'}, {'mean','std' }, 'DataVars', {'dff', 'tscore', 'respsnr'});

% get correlation within groups
dffcorr = splitapply(@(x,y) corr(x,y), double(stattab.respsnr), double(stattab.cellstim), findgroups(stattab.cell));





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



