function [I] = momentofinertia(input)
%MOMENTOFINERTIA Moment Of Inertia
%   Calculates Moment Of Inertia of given 2D or 3D object around its Center Of Mass
%   as the second moment of mass with respect to distance from a vertical
%   axis through the center of mass.
%
%   If input is 3D the MoI is calculated for each slice.
%
%   Input can be 2D or 3D matrix of reals, integers or logicals
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   01/06/2017
%   Last update:  11/04/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    %% get center of mass
    [cx, cy] = centerofmass(input);
    
    %% get Moment Of Inertia
    input = double(input);
    switch ndims(input)
        case 3
            % slicewise moment of inertia
            I = zeros(1, size(input,3));
            coor(:,1) = repmat([1:size(input,2)]', size(input,1), 1);
            coor(:,2) = reshape(repmat([1:size(input,1)], size(input,2), 1), numel(input(:,:,1)), 1);

            for slice = 1:size(input,3)
                C = [cx(slice) cy(slice)];
                dist = pdist2(coor, C)';
                mass_rows = reshape(input(:,:,slice)', numel(input(:,:,slice)), 1);
                I(slice) = (dist.^2)*mass_rows;
            end
        case 2
            % get coordinates array
            coor(:,1) = repmat([1:size(input,2)]', size(input,1), 1);
            coor(:,2) = reshape(repmat([1:size(input,1)], size(input,2), 1), numel(input), 1);

            C = [cx cy];
            % euclidean distance from the center of mass
            dist = pdist2(coor, C)';
            mass_rows = reshape(input', numel(input), 1);
            % MoI
            I = (dist.^2)*mass_rows;
    end
    
end