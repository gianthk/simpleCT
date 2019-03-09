function [BW1] = endosteummask(BW, closepixels, auto, dil_steps)
%ENDOSTEUMMASK get endosteum contour mask from binary mask of cortical bone
%   [BW1] = endosteummask(BW, closepixels, auto)
%   First of all an imclose-imopen step is applied in order to
%   remove all sparse bone pixels. Erosion and Dilation are performed
%   with 2D disk structuring element if the given mask is 2D,
%   sphere element if the given mask is 3D. Default radius = 3 pixels
%   ENDOSTEUMMASK attempts than to trace the boundary between bone and background starting
%   from the cortex midpoint (calculated as the mask center of mass).
%
%   [BW1] = endosteummask(BW) automatic endosteum
%   mask detection; initial erode-dilate step performed with
%   structuring element of radius = 3 pixels (Default)
%
%   [BW1] = endosteummask(BW, 5) automatic endosteum
%   mask detection; initial imclose-imopen step performed with
%   structuring element of radius = 5pixels
%
%   [BW1] = endosteummask(BW, [], 0) manual endosteum
%   mask drawing. Only available for 2-D.
%
%   [BW1] = endosteummask(BW, 'regionprops') endosteum mask detection with
%   regionprops method. The endosteum is considered as the pore surrounded
%   by bone tissue with larger area. Only available for 2-D.
%
%   Class Support
%   -------------
%   BW can be 2-D or 3-D matrix of logicals.
%
%   See also BWTRACEBOUNDARY, ENDOSTEUMCONTOUR
%   ______________________________________________________
%
%   Author: Gianluca Iori       <gianthk.iori@gmail.com>
%   BSRT - Charite Berlin
%   Created on:   01/11/2016
%   Last update:  09/03/2019
%
%   this function is part of the synchro toolbox
%   ______________________________________________________

if nargin<4                 dil_steps = 1;                          end
if isempty(dil_steps)       dil_steps = 1;                          end
if nargin < 3,              auto = true;                            end
if nargin < 2,              closepixels = 3;                        end
if isempty(closepixels)     closepixels = 3;                        end

if isnumeric(closepixels)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUTOMATIC ENDOSTEUM DETECTION 

    % erode and dilate mask: remove sparse bone pixels
    if closepixels > 0
        se_voxels = round(closepixels/dil_steps);

        if size(BW,3)>1
            % Morphological structuring element: sphere, default radius = 3 pixels
            se3 = strel3d(2*se_voxels);
        else
            % Morphological structuring element: disk, default radius = 3 pixels
            se3 = strel('disk',se_voxels);
        end
        fprintf('closing mask.');
        for i=1:dil_steps
            BW = imdilate(BW,se3);
            fprintf('.');
        end
        for i=1:dil_steps
            BW = imerode(BW,se3);
            fprintf('.');
        end
        fprintf('\nopening mask.');
        for i=1:dil_steps
            BW = imerode(BW,se3);
            fprintf('.');
        end
        for i=1:dil_steps
            BW = imdilate(BW,se3);
            fprintf('.');
        end
        fprintf('\n');
    end

    % initialize output mask
    BW1 = ones(size(BW));

    if auto
        for i=1:size(BW1,3)
            % GET ENDOSTEUM CONTOUR
            contour = endosteumcontour(BW(:,:,i));

            %% create mask according to periosteum contour
            if ~isempty(contour)
                BW1(:,:,i) = ~roipoly(BW(:,:,i), contour(:,2), contour(:,1));
            end

        end
    else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MANUAL ENDOSTEUM DRAWING 
        if size(BW1,3) > 1
            error('Manual endosteum selection not implemented for 3D data!!');
        else
            % let the user draw manually a polygon of cortical-trabecular boundary and generates mask of trabecular core
            % plot data (resized)
            fprintf(['Rules for the segmentation of the cortical bone compartment according to Malo et al (Malo, M. K. H., et al. Longitudinal elastic properties and porosity of cortical bone tissue vary with age in human proximal femur. Bone 53.2 (2013): 451-458.)\n'...
                    '\t1) Cortical bone includes Haversian canals.\n'...
                    '\t2) If the size of a pore is less than twice the average size of pores at the region, then the pore is included in the cortex.\n'...
                    '\t3) For the larger pores, the endocortical boundary is set to split the pore.\n'...
                    '\t4) If the diameter of the pore is smaller than the distance from the pore to the endosteal region, then the pore is included in the cortex.\n']);
            figure('units','normalized','position',[0 0 1 1]);
            imagesc(BW); colormap gray; axis image;
            title('draw polyline at cortical-trabecular boundary! Double click last point to apply changes.');
            BW1 = ~roipoly;
            close;
        end
    end

    % convert to logical
    BW1 = logical(BW1);
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REGIONPROPS METHOD 
    
    method = validatestring(closepixels,{'regionprops'},mfilename,'CLOSEPIXELS(METHOD)',1);
    
    % initialize output image
    BW1 = true(size(BW));
    
    % remove the background
    BW = imfill(BW, [1,1]);
   
    % track the pores
    props = regionprops(~BW, 'Area', 'PixelIdxList');
    for i=1:length(props)
        A(i) = props(i).Area;
    end        
    
    % find the pore with max area
    [val, id] = find(A==max(A));
    BW1(props(id).PixelIdxList)=false;

end