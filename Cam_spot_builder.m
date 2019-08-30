function [ Xcoordinates , Ycoordinates ] = Cam_spot_builder( pattern, sizesave, xsave ,ysave  )
%UNTITLED2 Summary of this function goes here
%   This function builds the x and y coordinates of the target spots in the
%   camera plane (different from the SLM plane)

for idx1 = 1:size(xsave,2);
    xmat=(1:size(pattern,1))-xsave(idx1); ymat=(1:size(pattern,1))-ysave(idx1);
    [Xmat,Ymat]=meshgrid(xmat,ymat);
    [Xcoordinates{idx1},Ycoordinates{idx1}]=find((Xmat.^2+Ymat.^2)<=(sizesave(idx1)/2).^2);
    if sizesave(idx1)==1
        Xcoordinates{idx1}=ysave(idx1); Ycoordinates{idx1}=xsave(idx1);
    end
        
    Xcoordinates{idx1}=round(Xcoordinates{idx1}); Ycoordinates{idx1}=round(Ycoordinates{idx1});


end

