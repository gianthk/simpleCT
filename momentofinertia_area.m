function [Ixx, Iyy] = momentofinertia_area(input)
%MOMENTOFINERTIA_AREA Area Moment Of Inertia
%   Calculates Second Moment(s) Of Area (or Area Moment(s) Of Inertia) 
%   around x-axis (Ixx) and y-axis (Iyy) of a given 2D binary image describing an area surface.
%   X and y axes are considered to cross the input area through its center of mass.
%
%   The second moment of area Ixx is computed in Cartesian coordinates as:
%       Ixx = ∬ y^2 dxdy
%   similarly:
%       Iyy = ∬ x^2 dxdy 
%
%   If input is 3D the Area Moment Of Inertia is calculated for each slice.
%
%   Input can be 2D or 3D matrix of logicals
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   11/04/2018
%   Last update:  28/04/2021
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    validateattributes(input, {'logical'},{'nonempty'},mfilename,'input',1);

    %% get center of mass
    [cx, cy] = centerofmass(input);
    
    %% get Area Moment Of Inertia
    input = double(input);
    switch ndims(input)
        case 3
            % slicewise moment of inertia
            Ixx = zeros(1, size(input,3));
            Iyy = zeros(1, size(input,3));
            coor(:,1) = repmat([1:size(input,2)]', size(input,1), 1);
            coor(:,2) = reshape(repmat([1:size(input,1)], size(input,2), 1), numel(input(:,:,1)), 1);

            for slice = 1:size(input,3)
                % distance from the x-axis
                dist_x = coor(:,2) - cy(slice);

                % distance from the x-axis
                dist_y = coor(:,1) - cx(slice);

                % surface area to 1D
                mass_rows = reshape(input(:,:,slice), numel(input(:,:,slice)), 1);
    
                % moments of area
                Ixx(slice) = sum(mass_rows.*(dist_x.^2));
                Iyy(slice) = sum(mass_rows.*(dist_y.^2));
        
            end
            
        case 2
            % area coordinates array
            coor(:,2) = repmat([1:size(input,2)]', size(input,1), 1);
            coor(:,1) = reshape(repmat([1:size(input,1)], size(input,2), 1), numel(input), 1);
            
            % distance from the x-axis
            dist_x = coor(:,2) - cy;

            % distance from the x-axis
            dist_y = coor(:,1) - cx;

            % surface area to 1D
            mass_rows = reshape(input, numel(input), 1);
            
            % moments of area
            Ixx = sum(mass_rows.*(dist_x.^2));
            Iyy = sum(mass_rows.*(dist_y.^2));
            
    end
    
end
