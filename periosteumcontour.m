function [contour] = periosteumcontour(BW, startpoint)
% PERIOSTEUMCONTOUR coordinates of periosteum
%   contour = periosteumcontour(BW) gets (outer) periosteum contour from binary mask.
%   
%	INPUTS:     BW              binary mask (2D or 3D)
%               x_startpoint    col position for the starting point of the tracebowndary procedure
%                               Default: col coordinate of mask center of mass
%
%   OUTPUT:     contour         nx2 or nx3 coordinates matrix (row, col, slice)
%
%   Class Support
%   -------------
%   BW should be a 2D or 3D matrix of type logical.
%   ______________________________________________________
%
%   Author: Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   01/11/2016
%   Last update:  23/05/2018
%
%   See also BWTRACEBOUNDARY
%
%   this function is part of the synchro toolbox
%   ______________________________________________________

    graphics = false;

    if nargin < 2
        startpoint = [0 0];
        startpointreset = true;
    end
    
    if isempty(startpoint)
        startpoint = [0 0];
        startpointreset = true;
    end

    switch ndims(BW)
        case 2    
            % FIND INITIAL POINT FOR CONTOUR DETECTION 
            % get first non-zero pixel from top of image on center line (works only for vertical line)
            start_pixel = [0 0];

            % x_startpoint = round(size(BW,2)/2);
            if startpoint(1) == 0
                startpoint = round(centerofmass(BW));
                ypixel_count = 1;
            else
                ypixel_count = startpoint(2);
            end

            while start_pixel(1) == 0
                if BW(ypixel_count, startpoint(1)) == 1
                    start_pixel = [ypixel_count, startpoint(1)];
                else ypixel_count = ypixel_count + 1;
                end
                if ypixel_count > size(BW,1)     % if end of image is reached without encountering bone
                    start_pixel = [1, 1];           % exit while -> this will produce an empty output mask

                    % try again one pixel to the right
                    % [x_startpoint, y_midpoint] = centerofmass(BW(:,:,i));
                    % ypixel_count = 1;
                    % x_startpoint = x_startpoint + 1;
                end
            end

            %% GET PERIOSTEUM CONTOUR
            % WARNING:  contour will be expressed as (rows, columns) (in pixels) of image
            %           and the starting point will appear at the end of the contour array one more time!
            start_pixel = double(start_pixel);
            contour = bwtraceboundary(BW(:,:), start_pixel, 'E');          % ,8,Inf,'counterclockwise'
            
        case 3
            contour = [];
            for i=1:size(BW,3)
                % clear contour_s
                % FIND INITIAL POINT FOR CONTOUR DETECTION 
                % get first non-zero pixel from top of image on center line (works only for vertical line)
                start_pixel = [0 0];

                % x_startpoint = round(size(BW,2)/2);
                if startpoint(1) == 0
                    startpoint = round(centerofmass(BW(:,:,i)));
                    ypixel_count = 1;
                else
                    ypixel_count = startpoint(2);
                end

                while start_pixel(1) == 0
                    if BW(ypixel_count, startpoint(1), i) == 1
                        start_pixel = [ypixel_count, startpoint(1)];
                    else ypixel_count = ypixel_count + 1;
                    end
                    if ypixel_count > size(BW(:,:,i),1)     % if end of image is reached without encountering bone
                        start_pixel = [1, 1];               % exit while -> this will produce an empty output mask
                    end
                end

                %% GET PERIOSTEUM CONTOUR
                % WARNING:  contour will be expressed as (rows, columns) (in pixels) of image
                %           and the starting point will appear at the end of the contour array one more time!
                start_pixel = double(start_pixel);
                contour_s = bwtraceboundary(BW(:,:,i), start_pixel, 'E');          % ,8,Inf,'counterclockwise'
                if ~isempty(contour_s)
                    contour = [contour; [contour_s ones(size(contour_s,1),1)*i]];              % append 2D contour points
                end
                
                if startpointreset
                    startpoint = [0 0];
                end
            end
            if graphics
                figure; plot3(contour(:,1), contour(:,2), contour(:,3));
            end
    end
end