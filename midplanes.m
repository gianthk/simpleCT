function [XY, XZ, YZ] = midplanes(I, varargin)
%MIDPLANES volume data midplanes images
%   [XY, XZ, YZ] = MIDPLANES(I) midplanes images of volume data I.
%
%   Example 1
%   -------
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   03/12/2017
%   Last update:  06/07/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

switch ndims(I)
    case 3
        slices = round(size(I)/2);
        
        XY = I(:,:,slices(3));
        YZ = squeeze(I(:,slices(2),:));     YZ = YZ';
        XZ = squeeze(I(slices(1),:,:));     XZ = XZ';
        
    case 2
        XY = I;
        XZ = [];
        YZ = [];
end

XY = double(XY);
XZ = double(XZ);
YZ = double(YZ);

end