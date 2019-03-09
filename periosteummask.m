function [BW1] = periosteummask(BW, closepixels, pad, dil_steps, a_thresh_pix, x_startpoint)
%PERIOSTEUMMASK calculates mask of periosteum
%   [BW1] = periosteummask(BW) gets periosteum contour mask from binary mask of cortical bone BW
%   [BW1] = periosteummask(BW, closepixels) control size of structuring element for dilate-erode steps
%   [BW1] = periosteummask(BW, [], pad) activate image zeros-padding (default = off)
%   [BW1] = periosteummask(BW, 10, [], 4) control number of dilate-erode sub steps
%   [BW1] = periosteummask(BW, 10, [], [], 4) control area threshold for pre-cleaning step of non connected struts
%   [BW1] = periosteummask(BW, [], [], [], [], x_startpoint) control the starting column for vertical search of the periosteum boundary
%   
%	INPUTS:
%       BW              @type logical       binary mask (2D or 3D)
%
%       closepixels     @type numeric       size of structuring element for dilate-erode step
%                                           (needed for filling of problematic holes)
%                                           Default: 3 pixels
%                                           If 0 is given no dilate-erode procedures is performed
%                     
%       x_startpoint    @type numeric       x position for the starting point of the tracebowndary procedure
%                                           Default: X coordinate of mask center of mass
%   ______________________________________________________
%
%   Author: Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   01/11/2016
%   Last update:  27/02/2018
%
%   See also BWTRACEBOUNDARY
%
%   this function is part of the synchro toolbox
%   ______________________________________________________

    if nargin<6                 x_startpoint = [];                      end
    if nargin<5                 a_thresh_pix = [];                      end
    if nargin<4                 dil_steps = 1;                          end
    if isempty(dil_steps)       dil_steps = 1;                          end
    if nargin<3                 pad = false;                            end
    if isempty(pad)             pad = false;                            end
    if nargin<2                 closepixels = 3;                        end
    if isempty(closepixels)     closepixels = 3;                        end

    %% ================ (0) CLEANING ================
    % remove sparse structs around mask
    if ~isempty(a_thresh_pix)

        for i = 1:size(BW,3)
            BWslice = BW(:,:,i);
            CC = bwconncomp(BWslice);
            numPixels = cellfun(@numel,CC.PixelIdxList);
            idx = find(numPixels < a_thresh_pix);

            for j = 1:length(idx)
                BWslice(CC.PixelIdxList{idx(j)}) = 0;
            end
            BW(:,:,i) = BWslice;

        end
        clear BWslice
    end
    
    %% add padding
    if pad
        BW = padarray(BW,[100 100 0]);
    end
    
    %% dilate and erode mask: fill eventual holes
    if closepixels > 0
        se_voxels = round(closepixels/dil_steps);
        
        if size(BW,3)>1
            % Morphological structuring element: sphere, default radius = 3 pixels
            % se3 = strel('sphere',closepixels);        % 'sphere' method not available in all matlab distributions!
            se3 = strel3d(se_voxels*2);
        else
            % Morphological structuring element: disk, default radius = 3 pixels
            se3 = strel('disk',se_voxels);
        end
        fprintf('dilating mask...\n');
        for i=1:dil_steps
            BW = imdilate(BW,se3);
        end
        
    end

    %% initialize output mask
    BW1 = ones(size(BW));

    for i=1:size(BW1,3)
        %% get periosteum contour
        contour = periosteumcontour(BW(:,:,i), x_startpoint);
        
        %% create mask according to periosteum contour
        if ~isempty(contour)
            BW1(:,:,i) = roipoly(BW(:,:,i), contour(:,2), contour(:,1));
        end

    end

    %% convert to logical
    BW1 = logical(BW1);

    %% erode back
    if closepixels > 0
        fprintf('eroding mask...\n');
        for i=1:dil_steps
            BW1 = imerode(BW1,se3);
        end
        
    end
    
    % remove padding
    if pad
        BW1 = BW1(101:end-100,101:end-100,:);
    end

end