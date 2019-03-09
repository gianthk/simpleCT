function BW1 = removenonconnectedstruts(BW, method, minarea)
% REMOVENONCONNECTEDSTRUTS remove all struts not connected to main one
    %   BW1 = removenonconnectedstruts(BW) uses regionprop3 (only available in Matlab > R2017b).
    %   INPUTS:     BW @type logical                binary mask (2D or 3D)
    %               SE @type numeric or strel       radius or structuring element (Default disk with radius = 50 pixels)
    %               METHOD @type char               'regionprops3', 'imclose' (Default), 'imfill', 'regionprops'
    %
    %   OUTPUT:     BW1 @type logical       output binary mask (2D or 3D)
	%
	%	METHODS:	'imclose' (Default)		while loop of imclose	(imdilate followed by imerode)
	%				'imfill'				imfill(BW, 'holes')		(slicewise if BW is 3D)
	%				'regionprops'			remove regions not belonging to main strut
    %				'regionprops2'			slicewise remove regions not belonging to main strut
    %				'area'                  remove regions based on minimum area (slicewise for 3D)
	%
    %   imclose method courtesy of Daniel Rohrbach
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   22/11/2018
    %   Last update:  28/01/2019
    %
    %   this function is part of the synchro toolbox    
    %   ______________________________________________________
    
    if nargin < 3,          minarea = 10;   end
    if isempty(minarea),  	minarea = 10;   end
    
    if nargin < 2,          method = 'imclose';     end
    if isempty(method),     method = 'imclose';     end
    
    method = validatestring(method,{'imclose','imfill','regionprops','regionprops2','area'},mfilename,'METHOD',1);
    
    switch method
        case 'imclose'
            warning('not implemented');
            return
            if isnumeric(se)
                se = strel('disk', se);
            end

            N0=length(find(BW));
            N1=N0+1;
            while N1 > N0
                N0=N1;
                BW1 = imclose(BW,se);
                N1=length(find(BW1));
            end
               
        case 'imfill'
            warning('not implemented');
            return
            
            BW1 = false(size(BW));
            for i = 1:size(BW,3)
                BW1(:,:,i) = imfill(BW(:,:,i), 'holes');
            end
            
        case 'regionprops'
            switch ndims(BW)
                case 3
                    if exist('regionprops3')==0
                        error('regionprops3 method is available only for Matlab > R2017b.. use regionprops2 instead');
                        return
                    end
                    BW1 = BW;
                    RP = regionprops3(BW1, 'Volume', 'VoxelIdxList');

                    boneV = max(RP.Volume);
                    for i=1:length(RP.Volume)
                        if RP.Volume(i) ~= boneV
                            BW1(RP.VoxelIdxList{i}) = 0;
                        end
                    end
                case 2
                    BW1 = BW;
                    RP = regionprops(BW1, 'Area', 'PixelIdxList');
                    Area = zeros(length(RP), 1);
                    for i=1:length(RP)
                        Area(i) = RP(i).Area;
                    end

                    boneA = max(Area);
                    for i=1:length(RP)
                        if RP(i).Area ~= boneA
                            BW1(RP(i).PixelIdxList) = 0;
                        end
                    end
            end
            
        case 'regionprops2'
            BW1 = BW;
            for i = 1:size(BW,3)
                BWslice = BW(:,:,i);
                RP = regionprops(BWslice, 'Area', 'PixelIdxList');
                Area = zeros(length(RP), 1);
                for j=1:length(RP)
                    Area(j) = RP(j).Area;
                end
                boneA = max(Area);
                for j=1:length(RP)
                    if RP(j).Area ~= boneA
                        BWslice(RP(j).PixelIdxList) = 0;
                    end
                end
                BW1(:,:,i) = BWslice;
            end
            
        case 'area'
            BW1 = BW;
            for i = 1:size(BW,3)
                BWslice = BW(:,:,i);
                CC = bwconncomp(BWslice);
                numPixels = cellfun(@numel,CC.PixelIdxList);
                idx = find(numPixels < minarea);

                for j = 1:length(idx)
                    BWslice(CC.PixelIdxList{idx(j)}) = 0;
                end
                BW1(:,:,i) = BWslice;
            end
         
    end
end
