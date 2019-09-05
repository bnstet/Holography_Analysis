function  save_name = run_extract_from_expt_table(tab_file, base_dir, row_num)

%%% NOTE: base_dir should be at the base of the mouse_name/date/... tree!
%%% e.g. videos are first found in base_dir/mouse_name/date/

disp(sprintf("Pulling row %d from table %s in base directory %s", row_num, tab_file, base_dir));

if ~isfile(tab_file)
    error(strcat("Input file ", tab_file, " does not exist!"))
end

mytab = readtable(tab_file, 'Delimiter', ',');
mytab = mytab(row_num,:);
date = int32(mytab.date);
mouse = mytab.mouse{:};

disp(sprintf("Mouse %s : date %d", mouse, date));

mddir = fullfile(base_dir, mouse, string(date));
mov_dir = fullfile(mddir, 'aligned');

mov_pattern = tif_pattern_from_initial(mytab.tif{:});

disp(sprintf("Looking for movie file pattern %s in directory %s", mov_pattern, mov_dir));

holo_path = char(fullfile(mddir,mytab.hologram{:}));
mask_dir = char(mytab.maskdir);
if mask_dir==char(NaN)
    mask_dir = 'masks';
end
if ~isfolder(fullfile(mddir,mask_dir))
    error(strcat("Can't find masks file ", mask_dir, " in directory ", mddir))
end
mask_dir = char(fullfile(mddir,mask_dir));

expt_info_file = char(fullfile(mddir,mytab.voyeur{:}));

red_channel = 0;

% primary extraction loop
save_name = traces_from_files(mov_dir, mov_pattern, holo_path, mask_dir, expt_info_file,red_channel);


save(save_name, 'date', 'mouse', '-append');

end