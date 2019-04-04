classdef registration3D < handle
%   sct class for 3D registration of 2D slice within 3D data volume
%
%   this class implements methods for slicewise registration of all the
%   slices of a given 3D dataset with a 2D image (generally a section of the
%   3D volume) acquired with a different modality and resolution.

%   once the slices of the 3D stack have been registered with the 2D image,
%   the correlation3D method allows to find the sectioning plane on which the 2D slice is laying.
%   with MATLAB's imwarp or with the sct package data3D can than be rotated
%   and sectioned in order to extract from it the very same slice as data2D
%   ______________________________________________________
%
%   Author: Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   18/08/2015
%   Last update:  25/01/2018
%
%   this class is part of the synchro toolbox    
%   ______________________________________________________

properties
    data2D              % (sctdata)
                        % the FIXED image
                        % 2D matrix of slice pixeldata
                        % [rows columns] of the matrix correspond to
                        % [Y X] dimensions of the dataset

    data2Dscaled        % (sctdata)
                        % image data2D scaled to fit data3D slice size.
                        % [rows columns] of the matrix correspond to data
                        % [Y X] dimensions.
                        
    data3D              % (sctdata)
                        % the MOVING image
                        % 3D matrix of volume data
                        % [rows columns slices] of the matrix correspond to
                        % [Y X Z] dimensions of the dataset
                        
    data3Dscaled        % (sctdata)
                        % WARNING: scaling 3D may cause memory overflow
                        % 3D matrix of volume data scaled to fit data2D slice size.
                        % [rows columns slices] of the matrix correspond to
                        % [Y X Z] dimensions of the dataset.
                        
    data3Dreg           % (sctdata)
                        % 3D matrix of registered volume data
                        % [rows columns slices] of the matrix correspond to
                        % [Y X Z] dimensions of the dataset
                        
    center              % (double) 1x3 array
                        % [x y z] coordinates of center of rotation
                        
    center_pix          % (double) 1x3 array of integers
                        % [row col slice] pixel coordinates of center of rotation
                        
    n                   % (double) 4x1 array
                        % fit plane coefficients in the form
                        %         n(1)*X + n(2)*Y + n(3)*Z + n(4) = 0
                        
    t                   % The distance threshold between data point and the plane
                        % used to decide whether a point is an inlier or
                        % not. default = 0.5 [mm].
                        
    regtype             % Image capture modality describes how your images have been captured,
                        % specified as either 'monomodal' (captured on the same device)
                        % or 'multimodal' (captured on different devices).
                        
    MaximumIterations   % Maximum number of slicewise iterations for optimizer (Default = 20)
    
    InitialRadius       % Initial radius for optimizer (Default = 0.0065)
    
    optimizer           % (optimizer object)
                        % Optimization configuration describes the method for
                        % optimizing the similarity metric, returned as one of the optimizer objects,
                        % registration.optimizer.RegularStepGradientDescent or
                        % registration.optimizer.OnePlusOneEvolutionary
                        
    metric              % (metric object)
                        % Metric configuration describes the image similarity metric
                        % to be optimized during registration,
                        % returned as one of the metric objects,
                        % registration.metric.MeanSquares or
                        % registration.metric.MattesMutualInformation.
                        
    slicewise_tformtype % transformType: Geometric transformation to be applied to the image to be registered
                        % 'translation' | 'rigid' | 'similarity' | 'affine'
                        
    slicewise_tform     % Geometric transformation for slicewise registration,
                        % specified as a geometric transformation object, affine2d.
    
    reg2D_tform         % Geometric transformation for 2D registration,
                        % specified as a geometric transformation object, affine2d.
                        
    slicewise_R         % (double) 3x3
                        % 2D transformation matrix for the slicewise
                        % registration of data3D to data2D
                        
    regslice_n          % (int)
                        % slice number of the slicewise resistered stack
                        % with max(R2_pre)
                        
    R                   % (double) 3x3
                        % 3D transformation matrix that rotates the fit
                        % plane on which data2D is laying parallel to the
                        % x-y plane
                        
    regslice            % (sctdata)
                        % data3Dreg cross section corresponding to data2D image
                        
    corrcomp_image      % (double)
                        % collage image of local regions of data3D with maximum
                        % correlation with data2D
                        
    overlap_image       % (double)  
                        % RED-GREEN 2 layers image showing overlap of
                        % data2D and regslice
                        
    R2_pre              % (double) [size(data3D,3)x1]
                        % slicewise array of 2D correlation coefficients
                        % between data3D slice and data2D
                        
    R2_reg              % (double) [size(data3Dreg,3)x1]
                        % slicewise array of 2D correlation coefficients
                        % between data3Dreg slice and data2D
    
    xyz_maxcorrcoeff    % (double) 3xn array of [x y z] coordinates of regions 
                        % of maximum local correlation between data2D and data3D
                        
    n_corrROIs          % number of correlation3D ROIs. default = 10x10
    
    corroffset          % pixel offset for correlation 3D. [top bottom left right]
                        % perform correlation on a smaller portion of data. default = [0 0 0 0]
    
    graphics            % (logical) flag. true for graphical correlation and registration output
    
    interactive         % (logical) flag. true for interactive optimization of the slciewise registration
    
end

properties (Access = private)                    
    data3Dregscaled     % (sctdata)
                        % WARNING: scaling 3D may cause memory overflow
                        % 3D matrix of registered volume data scaled to fit data2D slice size
                        % [rows columns slices] of the matrix correspond to
                        % [Y X Z] dimensions of the dataset.
    
    rcs_maxcorrcoeff    % (double) 3xn array of [rows columns slices] coordinates of regions 
                        % of maximum local correlation between data2D and data3D.
                        
    P                   % The three points in the data set that were found to
                        % define a plane having the most number of inliers.
                        % The three columns of P defining the three points.
    
    inliers             % The indices of the points that were considered
                        % inliers to the fitted plane.

end

methods (Static = true, Access = private)
    
end

methods (Static = true)
    % slice-wise 2D correlation between data2D and each slice of a data3D
    function [rcs_maxcorrcoeff, corrcoeff3D, uCTim_composed, uCT3D_composed] = correlation3D(im, data3D, n_ROIs, pixeloffset, graphics)
    %[rcs_maxcorrcoeff, corrcoeff3D, uCTim_composed, uCT3D_composed] = correlation3D(im, data3D, n_ROIs, graphics)
    %   CORRELATION3D slice-wise 2D correlation between 2D image and each slice of a 3D dataset
    %
    %   [rcs_maxcorrcoeff] = correlation3D(slice, data3D, n_ROIs)
    %   Returns (n_ROIs x n_ROIs)x3 array of (row, col, slice) indexes of
    %   local regions of maximum correlation between 2D im and data3D
    %
    %   use optional argument pixeloffset [top bottom left right] to
    %   perform correlation on a smaller portion of data
    %
    %   the first two dimensions of data3D must fit the size of slice.
    %   Images are divided into n_ROIs x n_ROIs subregions.
    %   2D correlation (corr2) coeffitients are computed for the corresponding regions.
    %   graphics (logical) = true activates graphical output.
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   01/12/2015
    %   Last update:  14/03/2017
    %   ______________________________________________________

    if nargin < 5       graphics = false;           end                 % Default graphics = 'off'
    
    if nargin < 4       pixeloffset = [0 0 0 0];    end                 % Default pixeloffset = 0

    if nargin < 3       warning('n_ROIs not specified: 10 x 10 sectioning will be applied.');
                        n_ROIs = 10;                end                 % 10 x 10

    if isempty(pixeloffset)       pixeloffset = [0 0 0 0];    end       % Default pixeloffset = 0            
    
    uCTim_composed = [];
    uCT3D_composed = [];
    
    % scale in order to remove neg values
    im = image2uint8(im);

    % ROI rows and columns indexes arrays
    tmp = [1+pixeloffset(1):ceil((size(im,1)-pixeloffset(1)-pixeloffset(2))/(n_ROIs+1)):size(im,1)-pixeloffset(2)];
    rowind(:,1) = tmp(1:end-1)';                                % ROI start row
    rowind(:,2) = tmp(2:end)';                                  % ROI end row
    clear tmp;
    tmp = [1+pixeloffset(3):ceil((size(im,2)-pixeloffset(3)-pixeloffset(4))/(n_ROIs+1)):size(im,2)-pixeloffset(4)];
    colind(:,1) = tmp(1:end-1)';                                % ROI start col
    colind(:,2) = tmp(2:end)';                                  % ROI end col
    clear tmp;
    
    % ROI rows and columns centroids 
    c_ROIs = mean(colind,2);
    r_ROIs = mean(rowind,2);

    % initialize corr coefficient 3D space
    corrcoeff3D = zeros(n_ROIs, n_ROIs, size(data3D,3));        % initialize corr coefficient 3D space
    
    tic();
    disp('registration3D: starting correlation procedure..');
    fprintf('status:     0 %%');
            
    % for each slice of data3D: convert to int8
    nslices = size(data3D,3);
    for slice = 1:nslices
        % scale to uint8 in order to remove neg values
        uCTim = image2uint8(data3D(:,:,slice));

        % get 3D space of local correlation coefficient between slice and im
        for i=1:n_ROIs
            for j=1:n_ROIs
                corrcoeff3D(i,j,slice) = corr2(im(rowind(i,1):rowind(i,2),colind(j,1):colind(j,2)),uCTim(rowind(i,1):rowind(i,2),colind(j,1):colind(j,2)));
            end
        end
        stat = round(100*slice/nslices);
        fprintf(1,[repmat('\b',1,numel(num2str(stat))+1) '\b%i %%'], stat);
    end
    fprintf('\ncorrelation completed in %s sec.\n',num2str(toc));
    fprintf('finding regions of max correlation...\n');
    
    % calculate level with max corrcoeff3D for each ROI
    z_maxcorrcoeff = nan(n_ROIs);                               % initialize corrcoeff matrix with nans
    options = fitoptions('gauss2');                             % get fit options for gauss2 model
    options.Lower=[0 0 0 0 0 0];                                % modify lower boundary for gauss2 fit model (avoid negative valued fits)
    for i=1:n_ROIs
        for j=1:n_ROIs
            if sum(isnan(squeeze(corrcoeff3D(i,j,:))))==0,
                % fit correlation coefficient over z with gauss2 distribution
                f = fit([1:size(corrcoeff3D,3)]',permute(corrcoeff3D(i,j,:),[3 1 2]),'gauss2',options);
    %             plot(f,[1:size(corrcoeff3D,3)]',permute(corrcoeff3D(i,j,:),[3 1 2]));           % visualize fit results
                if (f.b1 < size(data3D,3) && f.b1 >= 1)         % assign only plausible z_levels of max corrcoeff
                    z_maxcorrcoeff(i,j) = round(f.b1);
                end
            end
        end
    end

    if graphics         %__________________________________________________
        % compose image from regions with max corrcoeff
%         z_subregion = [];
%         z_subregion_count = 1;
        uCTim_composed = int16(nan([rowind(end,2)-rowind(1,1) colind(end,2)-colind(1,1) 1]));
        minZ = min(min(z_maxcorrcoeff));
        uCT3D_composed = int16(nan([rowind(end,2)-rowind(1,1) colind(end,2)-colind(1,1) max(max(z_maxcorrcoeff))-minZ]));
        for i=1:n_ROIs
            for j=1:n_ROIs
                if isnan(z_maxcorrcoeff(i,j))
                    uCTim_composed(rowind(i,1)+1-rowind(1,1):rowind(i,2)-rowind(1,1),colind(j,1)+1-colind(1,1):colind(j,2)-colind(1,1)) = nan;
                else
                    uCTim_composed(rowind(i,1)+1-rowind(1,1):rowind(i,2)-rowind(1,1),colind(j,1)+1-colind(1,1):colind(j,2)-colind(1,1)) = data3D(rowind(i,1):rowind(i,2)-1,colind(j,1):colind(j,2)-1,z_maxcorrcoeff(i,j));
                    uCT3D_composed(rowind(i,1)+1-rowind(1,1):rowind(i,2)-rowind(1,1),colind(j,1)+1-colind(1,1):colind(j,2)-colind(1,1),z_maxcorrcoeff(i,j)-minZ+1) = data3D(rowind(i,1):rowind(i,2)-1,colind(j,1):colind(j,2)-1,z_maxcorrcoeff(i,j));
%                     z_subregion(z_subregion_count) = z_maxcorrcoeff(i,j);
%                     z_subregion_count = z_subregion_count + 1;
                end
            end
        end
        
%         uCT3D_composed = nan([size(uCTim_composed) 100]);
%         z_uCT3D_composed = min(z_subregion):((max(z_subregion)-min(z_subregion))/99):max(z_subregion);
        
        % plot it!
        figure;
        subplot(1,2,1);
        imagesc(uCTim_composed); colormap gray
        M = size(uCTim_composed,1);                     N = size(uCTim_composed,2);
        hold on;
        % plot ROI borders
        for k = 1:n_ROIs-1
            x = [1 N];                                  y = [rowind(k,2) rowind(k,2)]-rowind(1,1);
            plot(x,y,'Color','w','LineStyle','-');      plot(x,y,'Color','k','LineStyle',':');
            x = [colind(k,2) colind(k,2)]-colind(1,1);	y = [1 N];
            plot(x,y,'Color','w','LineStyle','-');      plot(x,y,'Color','k','LineStyle',':');
        end
        % plot slice number on top of it
        for i=1:n_ROIs
            for j=1:n_ROIs
                text(mean(colind(j,:)-colind(1,1)),mean(rowind(i,:))-rowind(1,1),num2str(z_maxcorrcoeff(i,j)),'Color','red');
            end
        end
        hold off
        
        subplot(1,2,2);
        imagesc(im(rowind(1,1):rowind(end,2), colind(1,1):colind(end,2))); colormap gray
        M = rowind(end,2);                              N = colind(end,2);
        hold on;
        for k = 1:n_ROIs-1
            x = [1 N];                                  y = [rowind(k,2) rowind(k,2)]-rowind(1,1);
            plot(x,y,'Color','w','LineStyle','-');      plot(x,y,'Color','k','LineStyle',':');
            x = [colind(k,2) colind(k,2)]-colind(1,1);	y = [1 N];
            plot(x,y,'Color','w','LineStyle','-');      plot(x,y,'Color','k','LineStyle',':');
        end
        hold off
    end                 %__________________________________________________

    % build rcs_maxcorrcoeff output
    rcs_maxcorrcoeff        = NaN(3, n_ROIs*n_ROIs);
    rcs_maxcorrcoeff(1,:)   = (repmat(r_ROIs,n_ROIs, 1))';
    rcs_maxcorrcoeff(2,:)   = reshape(repmat(c_ROIs,1, n_ROIs)', n_ROIs*n_ROIs, 1);
    rcs_maxcorrcoeff(3,:)   = reshape(z_maxcorrcoeff, n_ROIs*n_ROIs, 1);
    rcs_maxcorrcoeff(:,isnan(rcs_maxcorrcoeff(3,:))) = [];      % remove nans

    end
    
    % fit size (number of rows and cols) of sctdata1in to the size of sctdata2in
    function [sctdata1outdata] = fitsctdatasize(sctdata1in, sctdata2in)
        % function [sctdata1out] = fitsctdatasize(sctdata1in, sctdata1in)
        %       fit size (number of rows and cols) of sctdata1in to the size of sctdata2in
        %         
        %       if sctdata1in.data is 150x80x200 and sctdata2in.data is 200x120 fitsctdatasize will
        %       give sctdata1outdata of size 200x120x200
        %
        %   INPUTS:
        %       sctdata1in          (sctdata)
        %       sctdata2in          (sctdata)
        %
        %   OUTPUT:
        %       sctdata1outdata     (2D or 3D matrix)
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianthk.iori@gmail.com)
        %   BSRT - Charite Berlin
        %   Created on:   03/01/2017
        %   Last update:  03/01/2017
        %   ______________________________________________________

        if nargin < 2
            error('not enough input arguments..');
        end

        switch numel(sctdata1in.size)
            case 3
                [sctdata1outdata] = registration3D.pad3D(sctdata1in, sctdata2in);
            case 2            
                [sctdata1outdata] = registration3D.pad2D(sctdata1in, sctdata2in);
        end
    end
    
    function [sctdata1outdata] = pad2D(sctdata1in, sctdata2in)
        % function [sctdata1out] = pad2D(sctdata1in, sctdata1in)
        %       pad method for fitsctdatasize and 2D data
        %
        %   INPUTS:
        %       sctdata1in          (sctdata)
        %       sctdata2in          (sctdata)
        %
        %   OUTPUT:
        %       sctdata1outdata     (2D matrix)
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianthk.iori@gmail.com)
        %   BSRT - Charite Berlin
        %   Created on:   22/02/2018
        %   Last update:  22/02/2018
        %   ______________________________________________________
        
        % resize rows
        if sctdata1in.size(2) > sctdata2in.size(2)
            sctdata1outdata = sctdata1in.data(1+ceil((sctdata1in.size(2)-sctdata2in.size(2))/2):end+1-ceil((sctdata1in.size(2)-sctdata2in.size(2))/2),:);
        elseif sctdata1in.size(2) < sctdata2in.size(2)
            sctdata1outdata = padarray(sctdata1in.data, [ceil((sctdata2in.size(2)-sctdata1in.size(2))/2) 0]);
        else
            sctdata1outdata = sctdata1in.data;
        end

        if size(sctdata1outdata,1) > sctdata2in.size(2)       % 1 pixel too large
            sctdata1outdata = sctdata1outdata(1:end-1,:);
        end
        if size(sctdata1outdata,1) < sctdata2in.size(2)       % 1 pixel is missing
            sctdata1outdata = padarray(sctdata1outdata, [1 0]);
        end

        % resize cols
        if sctdata1in.size(1) > sctdata2in.size(1)
            sctdata1outdata = sctdata1outdata(:,1+ceil((sctdata1in.size(1)-sctdata2in.size(1))/2):end+1-ceil((sctdata1in.size(1)-sctdata2in.size(1))/2));
        elseif sctdata1in.size(1) < sctdata2in.size(1)
            sctdata1outdata = padarray(sctdata1outdata, [0 ceil((sctdata2in.size(1)-sctdata1in.size(1))/2)]);
        else
            sctdata1outdata = sctdata1in.data;
        end

        if size(sctdata1outdata,2) > sctdata2in.size(1)       % 1 pixel too large
            sctdata1outdata = sctdata1outdata(:,1:end-1);
        end

        if size(sctdata1outdata,2) < sctdata2in.size(1)       % 1 pixel is missing
            sctdata1outdata = padarray(sctdata1outdata, [0 1]);
        end
        
    end
    
    function [sctdata1outdata] = pad3D(sctdata1in, sctdata2in)
        % function [sctdata1out] = pad3D(sctdata1in, sctdata1in)
        %       pad method for fitsctdatasize and 3D data
        %
        %   INPUTS:
        %       sctdata1in          (sctdata)
        %       sctdata2in          (sctdata)
        %
        %   OUTPUT:
        %       sctdata1outdata     (3D matrix)
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianthk.iori@gmail.com)
        %   BSRT - Charite Berlin
        %   Created on:   22/02/2018
        %   Last update:  22/02/2018
        %   ______________________________________________________
        
        % resize rows
        if sctdata1in.size(2) > sctdata2in.size(2)
            sctdata1outdata = sctdata1in.data(1+ceil((sctdata1in.size(2)-sctdata2in.size(2))/2):end+1-ceil((sctdata1in.size(2)-sctdata2in.size(2))/2),:,:);
        elseif sctdata1in.size(2) < sctdata2in.size(2)
            sctdata1outdata = padarray(sctdata1in.data, [ceil((sctdata2in.size(2)-sctdata1in.size(2))/2) 0 0]);
        else
            sctdata1outdata = sctdata1in.data;
        end

        if size(sctdata1outdata,1) > sctdata2in.size(2)       % 1 pixel too large
            sctdata1outdata = sctdata1outdata(1:end-1,:,:);
        end

        if size(sctdata1outdata,1) < sctdata2in.size(2)       % 1 pixel is missing
            sctdata1outdata = padarray(sctdata1outdata, [1 0 0]);
        end

        % resize cols
        if sctdata1in.size(1) > sctdata2in.size(1)
            sctdata1outdata = sctdata1outdata(:,1+ceil((sctdata1in.size(1)-sctdata2in.size(1))/2):end+1-ceil((sctdata1in.size(1)-sctdata2in.size(1))/2),:);
        elseif sctdata1in.size(1) < sctdata2in.size(1)
            sctdata1outdata = padarray(sctdata1outdata, [0 ceil((sctdata2in.size(1)-sctdata1in.size(1))/2) 0]);
        else
            sctdata1outdata = sctdata1in.data;
        end

        if size(sctdata1outdata,2) > sctdata2in.size(1)       % 1 pixel too large
            sctdata1outdata = sctdata1outdata(:,1:end-1,:);
        end

        if size(sctdata1outdata,2) < sctdata2in.size(1)       % 1 pixel is missing
            sctdata1outdata = padarray(sctdata1outdata, [0 1 0]);
        end

        % resize slices
%         if sctdata1in.size(3) > sctdata2in.size(3)
%             sctdata1outdata = sctdata1outdata(:,:,1+ceil((sctdata1in.size(3)-sctdata2in.size(3))/2):end+1-ceil((sctdata1in.size(3)-sctdata2in.size(3))/2));
%         elseif sctdata1in.size(3) < sctdata2in.size(3)
%             sctdata1outdata = padarray(sctdata1outdata, [0 0 ceil((sctdata2in.size(3)-sctdata1in.size(3))/2)]);
%         else
%             sctdata1outdata = sctdata1in.data;
%         end
% 
%         if size(sctdata1outdata,3) > sctdata2in.size(3)       % 1 pixel too large
%             sctdata1outdata = sctdata1outdata(:,:,1:end-1);
%         end
% 
%         if size(sctdata1outdata,3) < sctdata2in.size(3)       % 1 pixel is missing
%             sctdata1outdata = padarray(sctdata1outdata, [0 0 1]);
%         end
        
    end
    
    % SLICEWISER2 calculates slicewise 2D image correlation coefficient
    function [R2] = slicewiseR2(stack, image)
    %function [R2] = slicewiseR2(stack, slice)
    %   REGISTRATION3D.SLICEWISER2 calculates slicewise 2D image correlation coefficient
    %   between each slice of stack and given image.
    %   returns [size(stack,3) x 1] array of slicewise correlation coefficients
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   17/11/2016
    %   Last update:  22/12/2017
    %
    %   See also    CORR2
    %   ______________________________________________________

        % resize smaller image to fit size of larger
        if size(image,1) ~= size(stack,1) || size(image,2) ~= size(stack,2)
            if size(image,1)*size(image,2) > size(stack,1)*size(stack,2)
                stack = registration3D.imresize(stack, [size(image) size(stack,3)]);
            else
                image = imresize(image, [size(stack,1) size(stack,2)]);
            end
        end

        % scale data to uint8
        image = image2uint8(image);
        stack = image2uint8(stack);

        % initialize R2
        R2 = zeros(size(stack,3),1);

        % calculate 2-D correlation coefficient for each slice
        for slice = 1:size(stack,3)
            R2(slice) = corr2(image, stack(:,:,slice));
        end

    end
    
    % registration3D.imresize resizes 2D or 3D data with given factor.
    function [imres] = imresize(im, factor)
    %function [imres] = imresize(im, factor)
    %   REGISTRATION3D.IMRESIZE resizes 2D or 3D data with given factor.
    %
    %   factor          can be either a scalar or [row col slice] size of
    %                   data after resize.
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   17/11/2016
    %   Last update:  22/12/2017
    %
    %   See also    IMRESIZE
    %   ______________________________________________________

        if length(factor) == 3
            if ndims(im) ~= 3       error('registration3D.imresize: ndims(im) must be equal to length(factor)');	end
            if factor(3) == size(im,3)
                % imresize slicewise
                fprintf('registration3D.imresize: resizing 3D data slicewise..\n');
                for slice = 1:size(im,3)
                    % resize 2D data and mask
                    imres(:,:,slice) = imresize(im(:,:,slice), factor(1:2));
                end
            else
                % interp3
                error('registration3D.imresize: 3D interpolation not implemented yet');
            end
        elseif length(factor) == 2
            if ndims(im) == 3
                % imresize slicewise
                fprintf('registration3D.imresize: resizing 3D data slicewise..\n');
                for slice = 1:size(im,3)
                    % resize 2D data and mask
                    imres(:,:,slice) = imresize(im(:,:,slice), factor(1:2));
                end
            else
                imres = imresize(im, factor);
            end
        elseif length(factor) == 1
            if ndims(im) == 3
                % interp3
                error('registration3D.imresize: 3D interpolation not implemented yet');
            elseif ismatrix(im)
                fprintf('registration3D.imresize: imresize 2D data..\n');
                imres = imresize(im, factor);
            end
        end

    end
    
end

methods
    % resize smaller image to fit the size of larger one (not recommended)
    function samesize(this)
    % function samesize(this)
    % NOT RECOMMENDED: use fitsize instead
    % resize smaller image (between data2D and data3D) to fit the size of larger one
    % ..creates data2Dscaled and data3Dscaled!
    % NOTE: samesize doesn't guarantee that after resize the images will
    %       have the same pixelsize but only the same number of rows and cols!
    %       in order to fit both the pixelsize as well as the number of rows and
    %       columns of data2D and data3D you should use fitsize instead!
    %
    % if this.data2D is 200x120 and this.data3D is 150x80x200 samesize will
    % write this.data3Dscaled of size 200x120x200
    
        warning('The use of registration3D.samesize is not recommended! After samesize, data2D and data3D will have the same the same number of rows and cols but different pixelsize. In order to fit both the pixelsize as well as the number of rows and columns of data2D and data3D you should use fitsize instead! ');
        
        if this.data2D.size(1)*this.data2D.size(2) > this.data3D.size(1)*this.data3D.size(2)
            % 2D data is larger tha 3D data
            this.data3Dscaled = this.data3D.resize([this.data2D.size([2 1]) this.data3D.size(3)]);
            this.data2Dscaled = this.data2D;
        elseif this.data2D.size(1)*this.data2D.size(2) < this.data3D.size(1)*this.data3D.size(2) || this.data2D.size(1) ~= this.data3D.size(1)
            % 2D data is smaller or have the same size (but different
            % aspect ratio) than data 3D
            this.data3Dscaled = this.data3D;
            this.data2Dscaled = this.data2D.resize(this.data3D.size([2 1]));
        else
            % 2D and 3D data have already the same size
            this.data3Dscaled = this.data3D;
            this.data2Dscaled = this.data2D;
        end
    end
    
    % fit pixelsize and size of data3Dscaled and data2Dscaled
    % if either data3Dscaled or data2Dscaled is empty fitsize creates them
    % from data3D and data2D
    function fitsize(this)
    % function fitsize(this)
    % fit pixelsize and size (number of rows and cols) of data2D and data3D
    % ..creates data2Dscaled and data3Dscaled!
    %
    % if this.data2D is 200x120 and this.data3D is 150x80x200 samesize will
    % write this.data3Dscaled of size 200x120x200
        
        if isempty(this.data3Dscaled.data) || isempty(this.data2Dscaled.data)isempty(this.data3Dscaled.data) || isempty(this.data2Dscaled.data)
            if this.data3D.voxelsize(1) > this.data2D.voxelsize(1)
                % resize to fit voxelsize of data2D
                this.data3Dscaled.data = imresize(this.data3D.data,(this.data3D.voxelsize(1)/this.data2D.voxelsize(1)));
                this.data3Dscaled.voxelsize([1 2]) = this.data2D.voxelsize([1 2]);
                this.data3Dscaled.voxelsize(3) = this.data3D.voxelsize(3);

                % data2D.data is phisically copied to data2Dscaled.data
                % this consumes little memory and allows us to play with the
                % size of data2Dscaled without modifying the original data2D
                this.data2Dscaled.data = this.data2D.data;
                this.data2Dscaled.voxelsize = this.data2D.voxelsize;

            elseif this.data3D.voxelsize(1) < this.data2D.voxelsize(1)
                warning('data2D.voxelsize > data3D.voxelsize.. fitsize not tested yet!');

                % resize to fit voxelsize of data2D
                this.data2Dscaled.data = imresize(this.data2D.data,(this.data2D.voxelsize(1)/this.data3D.voxelsize(1)));
                this.data2Dscaled.voxelsize([1 2]) = this.data3D.voxelsize([1 2]);

                % data3Dscaled is linked to data3D
                this.data3Dscaled = this.data3D;

            else
                % datascaled is linked to data
                this.data2Dscaled = this.data2D;
                this.data3Dscaled = this.data3D;
                
            end
        end
        
        % resize to fit size of data3Dscaled
        this.data2Dscaled.data = registration3D.fitsctdatasize(this.data2Dscaled, this.data3Dscaled);

    end
    
    function fitsize3D(this)
    % function fitsize(this)
    % fit pixelsize and size (number of rows and cols) of data2D and data3D
    % ..creates data2Dscaled and data3Dscaled!
    %
    % if this.data2D is 200x120 and this.data3D is 150x80x200 samesize will
    % write this.data3Dscaled of size 200x120x200
        
        if isempty(this.data3Dscaled.data) || isempty(this.data2Dscaled.data)isempty(this.data3Dscaled.data) || isempty(this.data2Dscaled.data)
            if this.data3D.voxelsize(1) > this.data2D.voxelsize(1)
                % resize to fit voxelsize of data2D
                this.data3Dscaled.data = imresize(this.data3D.data,(this.data3D.voxelsize(1)/this.data2D.voxelsize(1)));
                this.data3Dscaled.voxelsize([1 2]) = this.data2D.voxelsize([1 2]);
                this.data3Dscaled.voxelsize(3) = this.data3D.voxelsize(3);

                % data2D.data is phisically copied to data2Dscaled.data
                % this consumes little memory and allows us to play with the
                % size of data2Dscaled without modifying the original data2D
                this.data2Dscaled.data = this.data2D.data;
                this.data2Dscaled.voxelsize = this.data2D.voxelsize;

            elseif this.data3D.voxelsize(1) < this.data2D.voxelsize(1)
                warning('data2D.voxelsize > data3D.voxelsize.. fitsize not tested yet!');

                % resize to fit voxelsize of data2D
                this.data2Dscaled.data = imresize(this.data2D.data,(this.data2D.voxelsize(1)/this.data3D.voxelsize(1)));
                this.data2Dscaled.voxelsize([1 2]) = this.data3D.voxelsize([1 2]);

                % data3Dscaled is linked to data3D
                this.data3Dscaled = this.data3D;

            else
                % datascaled is linked to data
                this.data2Dscaled = this.data2D;
                this.data3Dscaled = this.data3D;
                
            end
        end
        
        % resize to fit size of data3Dscaled
        this.data3Dscaled.data = registration3D.fitsctdatasize(this.data3Dscaled, this.data2Dscaled);

    end
    
    function launchregistration3D(this)
    % launch 3D registration
    
        % reset to empty registration output
        this.rcs_maxcorrcoeff = [];
        this.xyz_maxcorrcoeff = [];
        this.R = [];
        this.center = [];
        
        % fit data2D or data3D in order for these to have same size and voxelsize
        if isempty(this.data2Dscaled.data)
            this.fitsize;
        end
        
        % launch correlation of 2D image within 3D data using an offset of [200 180 460 60] (pixels) (top-bottom-left-right)
        [this.rcs_maxcorrcoeff, coeff3D, this.corrcomp_image] = this.correlation3D(this.data2Dscaled.data, this.data3Dscaled.data, this.n_corrROIs, this.corroffset, this.graphics);
        clear coeff3D
        
        % [rows columns slices] to [X Y Z] coordinates
        this.xyz_maxcorrcoeff(1,:) = this.rcs_maxcorrcoeff(2,:)*this.data3Dscaled.voxelsize(1);
        this.xyz_maxcorrcoeff(2,:) = this.rcs_maxcorrcoeff(1,:)*this.data3Dscaled.voxelsize(2);
        this.xyz_maxcorrcoeff(3,:) = this.rcs_maxcorrcoeff(3,:)*this.data3Dscaled.voxelsize(3);

        % fit plane of max correlation coefficients with RANSACFITPLANE
        [this.n, this.P, this.inliers] = ransacfitplane(this.xyz_maxcorrcoeff, this.t, this.graphics);
        if this.n(3)<0,  this.n=-this.n;   end         % invert versor if this points to negative z
        if this.graphics
            disp('______________________________');
            disp('MAX correlation regions fitted with plane using RANSACFITPLANE..');
            disp(sprintf('plane versor: %f %f %f',this.n(1),this.n(2),this.n(3)));
            disp('______________________________');
        end
        
        % get rotation matrix that rotates n onto z-versor
        this.R = vectors2rotation3Dmatrix(this.n(1:3),[0 0 1]);
        
        % determine center of rotation
        this.center = mean(this.xyz_maxcorrcoeff([1 2],:),2);           % x and y coordinates as center of image slice
        %         n(1)*X + n(2)*Y + n(3)*Z + n(4) = 0
        this.center(3) = (-this.n(1)*this.center(1)-this.n(2)*this.center(2)-this.n(4))/this.n(3);
        this.center = this.center';

        % display rotated max corr points
        if this.graphics                                      % space of max crosscorr coefficients
            figure; scatter3(this.xyz_maxcorrcoeff(1,:),this.xyz_maxcorrcoeff(2,:),this.xyz_maxcorrcoeff(3,:),'.');
            xlabel('X-axis [mm]');  ylabel('Y-axis [mm]');  zlabel('Z-axis [mm]');
            for i=1:length(this.xyz_maxcorrcoeff)
                xyz_maxcorrcoeff_rotated(:,i) = (this.R*(this.xyz_maxcorrcoeff(:,i)-this.center'))+this.center';      % rotate xyz_maxcorrcoeff
            end
            hold on
            scatter3(this.center(1),this.center(2),this.center(3),'*r');
            scatter3(xyz_maxcorrcoeff_rotated(1,:),xyz_maxcorrcoeff_rotated(2,:),xyz_maxcorrcoeff_rotated(3,:),'.m');
            scatter3(this.xyz_maxcorrcoeff(1,this.inliers),this.xyz_maxcorrcoeff(2,this.inliers),this.xyz_maxcorrcoeff(3,this.inliers),'og');
            hold off
        end
        
    end
    
    % slicewise registration initialization
    function init_reg_slicewise(this, regtype, MaximumIterations, InitialRadius, tformtype, graphics)
        
        if nargin < 6,                  graphics = this.graphics;                       end
        if nargin < 5,                  tformtype = this.slicewise_tformtype;           end
        if nargin < 4,                  InitialRadius = this.InitialRadius;             end
        if nargin < 3,                  MaximumIterations = this.MaximumIterations;     end
        if nargin < 2,                  regtype = this.regtype;                         end
        
        if isempty(tformtype)           tformtype = this.slicewise_tformtype;           end
        if isempty(InitialRadius)       InitialRadius = this.InitialRadius;             end
        if isempty(MaximumIterations)	MaximumIterations = this.MaximumIterations;     end
        if isempty(regtype)             regtype = this.regtype;                         end
        
        this.slicewise_tformtype = tformtype;
        this.InitialRadius = InitialRadius;
        this.MaximumIterations = MaximumIterations;
        this.regtype = regtype;
        
        % Configurations for intensity-based registration
        % multimodal images (SAM and MicroCT)
        [this.optimizer, this.metric] = imregconfig(this.regtype);
        
        this.optimizer.MaximumIterations = this.MaximumIterations;
        this.optimizer.InitialRadius = this.InitialRadius;
        
        if graphics
            fprintf('Settings for slicewise 2D registration:\n');
            fprintf('\tmultimodal images (SAM and MicroCT)\n');
            fprintf('\tType:\t%s\n', this.regtype);
            fprintf('\tMax_iter:\t%i\n', this.optimizer.MaximumIterations);
            fprintf('\tInit_radius:\t%i\n', this.optimizer.InitialRadius);
            fprintf('\ttransform type:\t%s\n', this.slicewise_tformtype);
        end
        
    end
    
    % interactive optimization of the 2D registration settings
    function optim_reg_slicewise(this, slice)
        
        if nargin < 2,	slice = round(this.data3Dscaled.size(3)/2);	end
                
        retry = true;
    
        while retry
            % try 1: Intensity-based image registration
            % Affine transformation consisting of translation, rotation, scale, and shear.
            movingRegisteredDefault = imregister(this.data3Dscaled.data(:,:,slice), this.data2Dscaled.data, this.slicewise_tformtype, this.optimizer, this.metric);

            % check result
            reg = figure; imshowpair(this.data2Dscaled.data, movingRegisteredDefault)

            % proceed ?
            prompt = {'Maximum n Iterations:','Initial Radius:'};
            dlg_title = 'Proceed with registration parameters:';
            num_lines = 1;
            defaultans = {num2str(this.optimizer.MaximumIterations), num2str(this.optimizer.InitialRadius)};
            answer = inputdlg(prompt, dlg_title, num_lines, defaultans);
            MaximumIterations   = str2num(answer{1});
            InitialRadius       = str2num(answer{2});

            if (MaximumIterations ~= this.optimizer.MaximumIterations || InitialRadius ~= this.optimizer.InitialRadius)
                this.optimizer.MaximumIterations = MaximumIterations;
                this.MaximumIterations = MaximumIterations;
                this.optimizer.InitialRadius = InitialRadius;
                this.InitialRadius = InitialRadius;
                close(reg)
            else
                retry = false;
            end
        end
        
        fprintf('New settings for slicewise 2D registration:\n');
        fprintf('\tmultimodal images (SAM and MicroCT)\n');
        fprintf('\tType:\t%s\n', this.regtype);
        fprintf('\tMax_iter:\t%i\n', this.optimizer.MaximumIterations);
        fprintf('\tInit_radius:\t%i\n', this.optimizer.InitialRadius);
        fprintf('\ttransform type:\t%s\n', this.slicewise_tformtype);
        
    end
    
    % 2D slicewise registration
    function register_slicewise(this)
        % class support: the output slicewise registered stack will be
        % int16 (implement further class support)
        
        if ~strcmp(this.slicewise_tformtype, 'affine')
            fprintf('slicewise_tformtype should be set to affine..\n');
            fprintf('current transform type:\t%s\n', this.slicewise_tformtype);
            m = input('Do you want to continue, Y/N [Y]:','s')
            if m=='N'
                return
            end
        end
        
        fprintf('Starting slicewise 2D registration with settings:\n');
        fprintf('\tMax_iter:\t%i\n', this.optimizer.MaximumIterations);
        fprintf('\tInit_radius:\t%i\n', this.optimizer.InitialRadius);
        fprintf('\ttransform type:\t%s\n', this.slicewise_tformtype);
        
        tic();
        fprintf('status:     0 %%');
        nslices = this.data3Dscaled.size(3);
        for slice = 1:nslices
            this.data3Dreg.data(:,:,slice) = int16(imregister(this.data3Dscaled.data(:,:,slice), this.data2Dscaled.data, this.slicewise_tformtype, this.optimizer, this.metric));
            stat = round(100*slice/nslices);
            fprintf(1,[repmat('\b',1,numel(num2str(stat))+1) '\b%i %%'], stat);
        end
        fprintf('\nslicewise registration completed in %s sec.\n',num2str(toc));
        
        this.data3Dreg.voxelsize = [this.data2Dscaled.voxelsize([1 2]) this.data3Dscaled.voxelsize(3)];
    end
    
    % 3D rotation (NOT RECOMMENDED.. USE imwrap INSTEAD!!)
    function rotatedata3D(this)
        warning('The use of rotate3D is not recommended! use MATLABS imwarp instead! ');
        
        %         this.data3Dreg = sctdata;
        %         this.data3Dreg.data = sct.rotate3D(this.data3D.data, this.R, this.center_pix);
        %         this.data3Dreg.voxelsize = this.data3D.voxelsize;
    end
    
    function calcR2_pre(this, data3D)
        % R2_pre: slicewise 2D correlation coefficient between data2Dscaled and data3D
        if nargin < 2
            if isempty(this.data3D)
                warning('data3Dreg is empty! Run slicewise registration first..');
            else
                fprintf('calculating R2_pre array for data3Dreg..\n');
                this.R2_pre = registration3D.slicewiseR2(this.data3D.data, this.data2Dscaled.data);
            end
        else
            fprintf('calculating R2_pre array for input 3D data..\n');
            this.R2_pre = registration3D.slicewiseR2(data3D, this.data2Dscaled.data);
        end
    end
    
    function calcR2_reg(this, data3Dreg)
        % R2_reg: slicewise 2D correlation coefficient between data2Dscaled and data3Dreg
        if nargin < 2
            if isempty(this.data3Dreg)
                warning('data3Dreg is empty! Run slicewise registration first..');
            else
                fprintf('calculating R2_reg array for data3Dreg..\n');
                this.R2_reg = registration3D.slicewiseR2(this.data3Dreg.data, this.data2Dscaled.data);
            end
        else
            fprintf('calculating R2_reg array for input 3D data..\n');
            this.R2_reg = registration3D.slicewiseR2(data3Dreg, this.data2Dscaled.data);
        end
    end

    function getOverlap_image(this)
        this.overlap_image = imfuse(this.data2D.data, imresize(this.regslice.data,size(this.data2D.data)),'falsecolor','Scaling','independent','ColorChannels',[1 2 0]);
    end
    
    function getRegslice(this)
        if isempty(this.data3Dreg.data)     error('you must first register and rotate data3D!');    end
        if isempty(this.R2_reg)             error('you must first calculate R2 coefficient for reg data!');    end
        
        % plot R2 before and after registration and rotation
        if this.graphics
            figure;     plot(this.R2_pre);
            hold on;    plot(this.R2_reg);
            legend('before','registered');
        end
        
        % get slicenumber of slice with maximum correlation coefficient
        % after 3D rotation
        regslice_n = find(this.R2_reg == max(this.R2_reg));
        if length(regslice_n) > 1   regslice_n = mean(regslice_n);  end
        this.regslice.data = this.data3Dreg.data(:,:,regslice_n);
        this.regslice.voxelsize = this.data3D.voxelsize;
        this.regslice.headerfile = this.data3D.headerfile;
        this.regslice.rawfile = this.data3D.rawfile;
        this.regslice.elementtype = this.data3D.elementtype;
        
    end
    
    %% registration3D getters                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % size (in bytes) of the registration properties
    function GetSize(this)
        props = properties(this);
        totSize = 0;
        for ii=1:length(props)
            currentProperty = getfield(this, char(props(ii)));
            s = whos('currentProperty'); totSize = totSize + s.bytes;
        end
        fprintf(1, '%d bytes\n', totSize);
    end
    
    % pixel coordinates of the center of rotation for data3D
    function centerpixout = getCenter_pix(this, data3D, center)
        if nargin < 2
            if isempty(this.center)
                fprintf('getCenter_pix: center of rotation is empty! ..run correlation3D first.\n')
                centerpixout = [];
            elseif isempty(this.data3D.voxelsize)
                fprintf('getCenter_pix: data3D.voxelsize is empty! ..cannot determine center of rotation in pixels.\n')
                centerpixout = [];
            else
                centerpixout(1) = round(this.center(2)/this.data3D.voxelsize(2));       % rows
                centerpixout(2) = round(this.center(1)/this.data3D.voxelsize(1));       % cols
                centerpixout(3) = round(this.center(3)/this.data3D.voxelsize(3));       % slices
            end
        elseif nargin < 3
            fprintf('getCenter_pix: calculating center of rotation (in pixels) for input data3D..\n')
            if isempty(this.center)
                fprintf('getCenter_pix: center of rotation is empty! ..run correlation3D first.\n')
                centerpixout = [];
            elseif isempty(data3D.voxelsize)
                fprintf('getCenter_pix: input data3D.voxelsize is empty! ..cannot determine center of rotation in pixels.\n')
                centerpixout = [];
            else
                centerpixout(1) = round(this.center(2)/data3D.voxelsize(2));       % rows
                centerpixout(2) = round(this.center(1)/data3D.voxelsize(1));       % cols
                centerpixout(3) = round(this.center(3)/data3D.voxelsize(3));       % slices
            end
            elseif nargin < 4
            fprintf('getCenter_pix: calculating center of rotation (in pixels) for input data3D and center (in mm)..\n')
            if isempty(center)
                fprintf('getCenter_pix: input center of rotation is empty! ..run correlation3D first.\n')
                centerpixout = [];
            elseif isempty(data3D.voxelsize)
                fprintf('getCenter_pix: input data3D.voxelsize is empty! ..cannot determine center of rotation in pixels.\n')
                centerpixout = [];
            else
                centerpixout(1) = round(center(2)/data3D.voxelsize(2));       % rows
                centerpixout(2) = round(center(1)/data3D.voxelsize(1));       % cols
                centerpixout(3) = round(center(3)/data3D.voxelsize(3));       % slices
            end
        end
    end
    
    % best match slice before 3D registration
    function [regslice_n_pre, regslice_n_reg]  = getRegslice_n(this)
        regslice_n_pre = [];
        regslice_n_reg = [];
        
        fprintf('registration3D.getRegslice_n: no smooth implemented.\n');
        if ~isempty(this.R2_pre)
            regslice_n_pre = find(this.R2_pre == max(this.R2_pre));
        end
        
        if ~isempty(this.R2_reg)
            regslice_n_reg = find(this.R2_reg == max(this.R2_reg));
        end
    end
    
    function [results] = getSlicewise_results(this)
        results.MaximumIterations = this.MaximumIterations;
        results.InitialRadius = this.InitialRadius;
        results.R2_pre = this.R2_pre;
        results.maxR2_pre = max(this.R2_pre);
        results.regslice_n = this.regslice_n;
        results.slicewise_tform = this.slicewise_tform;
        results.reg2D_tform = this.reg2D_tform;
        
    end
    
    %% registration3D setters                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function setRegslice_n(this, prop)
        this.regslice_n = prop;
    end
    
	%% registration3D class constructor         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function this = registration3D()
        this.data2D = sctdata;
        this.data2Dscaled = sctdata;
        this.data3D = sctdata;
        this.data3Dreg = sctdata;
        this.data3Dscaled = sctdata;
        this.regslice = sctdata;
        this.corrcomp_image = [];
        this.overlap_image = [];
        this.rcs_maxcorrcoeff = [];
        this.xyz_maxcorrcoeff = [];
        this.center = [];
        this.center_pix = [];
        this.regtype = 'multimodal';
        this.MaximumIterations = 20;
        this.InitialRadius = 0.0065;
        this.optimizer = [];
        this.metric = [];
        this.slicewise_tformtype = 'affine';
        this.slicewise_tform = [];
        this.reg2D_tform = [];
        this.slicewise_R = [];
        this.regslice_n = [];
        this.R = [];
        this.n = [];
        this.P = [];
        this.inliers = [];
        this.t = 0.5;
        this.R2_pre = [];
        this.R2_reg = [];
        this.n_corrROIs = 10;
        this.corroffset = [];
        this.graphics = false;
        this.interactive = true;
                        
    end
end

end