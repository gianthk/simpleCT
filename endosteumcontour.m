function [contour] = endosteumcontour(BW)
% ENDOSTEUMCONTOUR Endosteum Contour
%   Finds point cloud defining endosteum contour from binary mask of cortical bone
%   
%	INPUTS:
%       BW              @type logical       binary mask (2D or 3D)
%   
%   contour will be expressed as (rows, columns) (in pixels) of image
%   and the starting point will appear at the end of the contour array one more time!
%   ______________________________________________________
%
%   Author: Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   01/11/2016
%   Last update:  25/04/2018
%
%   See also BWTRACEBOUNDARY
%
%   this function is part of the synchro toolbox
%   ______________________________________________________

    graphics = false;
    
    switch ndims(BW)
        case 2    
            % FIND INITIAL POINT FOR CONTOUR DETECTION
            % get first non-zero pixel from baricenter of image on vertical midline
            % the underlying hypothesis is that the cortical bone mask is hollow
            % and that the baricenter will fall into the bone marrow cavity
            % additionally, the contour detection will start the search of contour
            % points towards the right (EAST) of the mask.. thus supposing that the
            % endosteum is apporximatively circular. for cross sections containing
            % trabecular bone regions this may not be the case!

            start_pixel = [0 0];
            [x_midpoint, ypixel_count] = centerofmass(BW);
            x_midpoint = round(x_midpoint);
            ypixel_count = round(ypixel_count);

            if isnan(x_midpoint)
                x_midpoint = round(size(BW,2)/2);
            end

            if isnan(ypixel_count)
                ypixel_count = 1;
            end

            while start_pixel(1) == 0
                if BW(ypixel_count, x_midpoint) == 1
                    start_pixel = [ypixel_count, x_midpoint];
                else ypixel_count = ypixel_count + 1;
                end
                if ypixel_count > size(BW,1)     % if end of image is reached without encountering bone
                    start_pixel = [1, 1];           % exit while -> this will produce an empty output mask

                    % try again one pixel to the right
                    % [x_midpoint, ypixel_count] = centerofmass(BW(:,:,i));
                    % x_midpoint = x_midpoint + 1;
                end
            end

            %% GET ENDOSTEUM CONTOUR
            % WARNING:  contour will be expressed as (rows, columns) (in pixels) of image
            %           and the starting point will appear at the end of the contour array one more time!
            start_pixel = double(start_pixel);
            contour = bwtraceboundary(BW, start_pixel, 'E');          % ,8,Inf,'counterclockwise'
        case 3
            contour = [];
            for i=1:size(BW,3)
                clear contour_s
                % FIND INITIAL POINT FOR CONTOUR DETECTION
                % get first non-zero pixel from baricenter of image on vertical midline
                % the underlying hypothesis is that the cortical bone mask is hollow
                % and that the baricenter will fall into the bone marrow cavity
                % additionally, the contour detection will start the search of contour
                % points towards the right (EAST) of the mask.. thus supposing that the
                % endosteum is apporximatively circular. for cross sections containing
                % trabecular bone regions this may not be the case!

                start_pixel = [0 0];
                [x_midpoint, ypixel_count] = centerofmass(BW(:,:,i));
                x_midpoint = round(x_midpoint);
                ypixel_count = round(ypixel_count);

                if isnan(x_midpoint)
                    x_midpoint = round(size(BW(:,:,i),2)/2);
                end

                if isnan(ypixel_count)
                    ypixel_count = 1;
                end

                while start_pixel(1) == 0
                    if BW(ypixel_count, x_midpoint, i) == 1
                        start_pixel = [ypixel_count, x_midpoint];
                    else ypixel_count = ypixel_count + 1;
                    end
                    if ypixel_count > size(BW(:,:,i),1)     % if end of image is reached without encountering bone
                        start_pixel = [1, 1];           % exit while -> this will produce an empty output mask

                        % try again one pixel to the right
                        % [x_midpoint, ypixel_count] = centerofmass(BW(:,:,i));
                        % x_midpoint = x_midpoint + 1;
                    end
                end

                %% GET ENDOSTEUM CONTOUR
                % WARNING:  contour will be expressed as (rows, columns) (in pixels) of image
                %           and the starting point will appear at the end of the contour array one more time!
                start_pixel = double(start_pixel);
                contour_s = bwtraceboundary(BW(:,:,i), start_pixel, 'E');           % ,8,Inf,'counterclockwise'
                if ~isempty(contour_s)
                    contour = [contour; [contour_s ones(size(contour_s,1),1)*i]];              % append 2D contour points
                end
            end
            
            if graphics
                figure; plot3(contour(:,1), contour(:,2), contour(:,3));
            end
    end
end