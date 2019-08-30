function [Spotidx, SpotMat, spot_num, multistim] = load_spot_info(pattern_path)

% load spot info for 2p experiments; can work on directory of individual
% mat files or a single mat file
disp(sprintf("Loading stim spot info from %s", pattern_path));

if isfolder(pattern_path)
    multistim = true;
    fnamesPat = dir(fullfile(pattern_path,'*.mat'));
    % Check number of holograms, set up variables
    spot_num=length(fnamesPat);
    Spotidx{spot_num}=[];
    SpotMat=cell(spot_num,1);
    SpotMat(:)={zeros(512)};
    spots.xcoordsAll = [];
    spots.ycoordsAll = [];
    spots.sizesAll = [];

    % Cycle through each hologram and extract the coordinates. Turn these into
    % indices for linearly defining stimulated pixels in images.
    for idx1 = 1:spot_num

        % Load each pattern in directory
        fnamePat = fnamesPat(idx1).name;
        fpathPat = fnamesPat(idx1).folder;
        load(fullfile(fpathPat,fnamePat));% loading the pattern file

        spots.xcoordsAll(idx1) = xsave;
        spots.ycoordsAll(idx1) = ysave;
        spots.sizesAll(idx1) = sizesave;

        % Get coordinates for each pattern. Coordinates
        [Xcoordinates, Ycoordinates] = ...
            Cam_spot_builder(pattern, sizesave, xsave ,ysave);

        spots.XcoordinatesAll{idx1} = Xcoordinates;
        spots.YcoordinatesAll{idx1} = Ycoordinates;

        % Test for 1 pixel spots
        [c1,d1]=cellfun(@size ,Xcoordinates);
        [c2,d2]=cellfun(@size ,Ycoordinates);

        % if it is 1 pixels spot draws a circle around it for analysis
        if sum(c1.*d1.*c2.*d2)/size(c1,2)==1
            [ Xcoordinates ,Ycoordinates ] = ...
                CircleDrawer( Xcoordinates, Ycoordinates );
        end



        % Create lookup for spots
        Spotidx{idx1}=sub2ind(size(pattern),Xcoordinates{:},...
            Ycoordinates{:});
        SpotMat{idx1}(Spotidx{idx1})=1;


    end

else
    [~, ~, file_ext] = fileparts(pattern_path);
    multistim = false;
    if strcmp(file_ext,'.mat')
        load(pattern_path, 'pattern', 'sizesave', 'xsave', 'ysave');% loading the pattern file
        % Seperating the pattern to the different cells and show figs on subplots
        [ Xcoordinates , Ycoordinates ] = Cam_spot_builder( pattern, sizesave, xsave ,ysave  );
        spot_num=size(Xcoordinates,2);
        [c1,d1]=cellfun(@size ,Xcoordinates);[c2,d2]=cellfun(@size ,Ycoordinates);
        if sum(c1.*d1.*c2.*d2)/size(c1,2)==1  % if it is 1 pixels spot draws a circle around it for analysis
            [ Xcoordinates ,Ycoordinates ] = CircleDrawer( Xcoordinates, Ycoordinates );
        end
        Spotidx{spot_num}=[];
        SpotMat=cell(spot_num,1); SpotMat(:)={zeros(size(pattern))};
        for idx=1:spot_num
            Spotidx{idx}=sub2ind(size(pattern),Xcoordinates{idx},Ycoordinates{idx});
            SpotMat{idx}(Spotidx{idx})=1;
        end
    else
       Spotidx = cell(1,1); SpotMat = cell(1,1);
       temp = imread(pattern_path);
       [height, width] = size(temp);
       SpotMat{1} = logical(temp./max(max(temp)));
       Spotidx{1} = find(temp);
       spot_num = 1;
    end

end



end

