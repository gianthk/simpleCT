function BW1 = fillpores(BW, se, method)
% FILLPORES fills all pores with size smaller than given radius
    %   BW1 = fillpores(BW, se) applies while loop of imclose (imdilate followed by imerode)
    %   on input binary mask BW with structuring element se. ('imclose' METHOD)
	%	The loop stops when the mask remains unvaried before and after imclose step.
    %   INPUTS:     BW @type logical                binary mask (2D or 3D)
    %               SE @type numeric or strel       radius or structuring element (Default disk with radius = 50 pixels)
    %               METHOD @type char               'imclose' (Default), 'imfill', 'regionprops', 'regionprops3'
    %
    %   OUTPUT:     BW1 @type logical       output binary mask (2D or 3D)
    %
	%
	%	METHODS:	'imclose' (Default)		while loop of imclose	(imdilate followed by imerode)
	%				'imfill'				imfill(BW, 'holes')		(slicewise if BW is 3D)
	%				'regionprops'			close regions with MinorAxisLength < SE		(slicewise)
	%				'regionprops3'			close regions not belonging to main strut	(3D)
	%
    %   imclose method courtesy of Daniel Rohrbach
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   --/--/----
    %   Last update:  22/11/2018
    %
    %   this function is part of the synchro toolbox    
    %   ______________________________________________________
    
    if nargin < 3,          method = 'imclose';     end
    if isempty(method),     method = 'imclose';     end
    
    method = validatestring(method,{'imclose','imfill','regionprops','regionprops3'},mfilename,'METHOD',1);
    
    if nargin < 2,          se = strel('disk',50);  end
    if isempty(se)
        switch method
            case 'imclose'
                se = strel('disk',50);
            case 'regionprops'
                se = 50;
        end
    end
    
    switch method
        case 'imclose'
            if isnumeric(se)
                se = strel('disk', se);
            end

            N0=length(find(BW));
            N1=N0+1;
            while N1 > N0,
                N0=N1;
                BW1 = imclose(BW,se);
                N1=length(find(BW1));
            end
            
        case 'regionprops'
            BW1 = false(size(BW));
            for i = 1:size(BW,3)
                BWslice = ~BW(:,:,i);
                CC = bwconncomp(BWslice);
                RP = regionprops(BWslice, 'MinorAxisLength');
                th = cat(1, RP.MinorAxisLength);
                idx = find(th<2*se);

                for j = 1:length(idx)
                    BWslice(CC.PixelIdxList{idx(j)}) = 0;
                end
                BW1(:,:,i) = ~BWslice;
            end
            
        case 'regionprops3'
            if exist('regionprops3')==0
                error('regionprops3 method is available only for Matlab > R2017b.. use regionprops (slicewise) instead!');
                return
            end
            BW1 = ~BW;
            RP = regionprops3(BW1, 'Volume', 'VoxelIdxList');

            BGV = max(RP.Volume);
            for i=1:length(RP.Volume)
                if RP.Volume(i) < 0.1*BGV;
                    BW1(RP.VoxelIdxList{i}) = 0;
                end
            end

            BW1 = ~BW1;
            
        case 'imfill'
            BW1 = false(size(BW));
            for i = 1:size(BW,3)
                BW1(:,:,i) = imfill(BW(:,:,i), 'holes');
            end
    end
end
