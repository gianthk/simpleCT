function [cmassx_array, cmassy_array] = centerofmass(input, method, fill)
%CENTEROFMASS Center of Mass
%   [cmassx, cmassy] = centerofmass(input) calculates baricenter coordinates
%   using x and y projection method
%
%   [cmassx, cmassy] = centerofmass(input, 'mass') calculates baricenter
%   coordinates using center of mass method
%
%   [cmassx, cmassy] = centerofmass(input, '', 1) calculates baricenter
%   coordinates using x and y projection method. the (binary) input image
%   is filled first
%
%   Baricenter coordinates cmassx_array, cmassy_array are returned in pixels.
%   If input is 3D matrix centerofmass treats the 3rd dimension as Zaxis
%   and calculates baricenter slicewise, returning arrays of baricenter
%   coordinates both of length = size(input,3).
%
%   Input can be 2D or 3D matrix of reals, integers or logicals
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   30/08/2016
%   Last update:  06/04/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    if(~exist('method', 'var')),    method = 'mass';            end
    if(isempty(method)),            method = 'mass';            end
    
    method = validatestring(method, {'projection','mass'},mfilename,'METHOD',2);

    if(~exist('fill', 'var')),      fill = 0;                   end

    %% fill input binary image
    if(fill)
        if(~islogical(input))
            warning('imfill skipped: input is not logical!');
        else
            input = imfill(bwareaopen(input,10000),'holes');
        end
    end

    %% get center of mass
    switch islogical(input)
        case true
            switch method
                case 'mass'
                    for kk=1:size(input,3)
                        cmassx = 0;
                        cmassy = 0;
                        pixels = 0;

                        % y-axis should correspond to the first dimension of the matrix
                        % and x-axis should correspond to the second dimension of the
                        % matrix
                        for ii=1:size(input,1)
                            for jj=1:size(input,2)
                                if(input(ii,jj,kk) == 1)
                                    cmassy = cmassy + ii;
                                    cmassx = cmassx + jj;
                                    pixels = pixels + 1;
                                end
                            end
                        end
                        cmassx_array(kk) = cmassx / pixels;
                        cmassy_array(kk) = cmassy / pixels;
                    end

                case 'projection'
                    for kk=1:size(input,3)
                        y = sum(input(:,:,kk),2);
                        cmassy = sum(y'*[1:length(y)]');
                        cmassy_array(kk) = cmassy/sum(y);

                        x = sum(input(:,:,kk),1);
                        cmassx = sum(x*[1:length(x)]');
                        cmassx_array(kk) = cmassx/sum(x);
                    end
            end
            
        case false
            input = double(input);
            switch method
                case 'mass'
                    x_coor = repmat([1:size(input,2)]', size(input,1), 1);
                    y_coor = repmat([1:size(input,1)]', size(input,2), 1);
                    
                    for kk=1:size(input,3)
                        mass_rows = reshape(input(:,:,kk), 1, numel(input(:,:,kk)));
                        mass_cols = reshape(input(:,:,kk)', 1, numel(input(:,:,kk)));
                        cmassy_array(kk) = (mass_rows*y_coor) / sum(mass_rows);
                        cmassx_array(kk) = (mass_cols*x_coor) / sum(mass_cols);
                    end

                case 'projection'
                    % error('projection method requires input data to be binary!');
                    for kk=1:size(input,3)
                        y = sum(input(:,:,kk),2);
                        cmassy = sum(y'*[1:length(y)]');
                        cmassy_array(kk) = cmassy/sum(y);

                        x = sum(input(:,:,kk),1);
                        cmassx = sum(x*[1:length(x)]');
                        cmassx_array(kk) = cmassx/sum(x);
                    end
            end

end