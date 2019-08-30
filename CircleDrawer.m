function [ Xcoordinates ,Ycoordinates ] = CircleDrawer( Xcoordinates, Ycoordinates )
%UNTITLED4 Summary of this function goes here
%   This function draws a circle around 1 pixels spot in order to analyze
%   them
radCurve=str2num(cell2mat(inputdlg('Which radius [pixels] around the points do you want to analyze ?','Analysis radius')));
imSize = 512;

for idx=1:size(Xcoordinates,2)
    
    x = [];
    y = [];
    
    for ii=-imSize/2:imSize/2
        for jj = -imSize/2:imSize/2
            if ( sqrt(ii^2 + jj^2) < radCurve )
                x = [x ii];
                y = [y jj];
            end
        end
    end
    
    Xcoordinates{idx} = floor(x)+Xcoordinates{idx};
    Ycoordinates{idx} = floor(y)+Ycoordinates{idx};
end

end

