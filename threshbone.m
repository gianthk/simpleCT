function [BW, level, method] = threshbone(I, varargin)
%THRESHBONE synchro thresholding function for bone images
%   BW = THRESHBONE(I, METHOD) thresholds input image I with given METHOD.
%   Possible values for METHOD are:
%
%     'otsu'        (Default) Otsu's method.
%     'fixed'       Single or double level fixed global threshold.
%     'lakshmanan'  Adaptive threshold. (use for SAM images) (Lakshmanan_2007)
%     'ridler'      Global threshold using iterative isodata method. (Ridler_1987)
%     'adaptive'    Local Adaptive Threshold (LAD) mean or median based. (Wellner_1993)
%
%   Depending on METHOD, THRESHBONE may take one additional parameter
%   which you can supply.
%
%   BW = THRESHBONE(I, 'otsu') thresholds input image I with (single level) Otsu's method.
%   BW is a mask of logicals with the same size as I.
%
%   [BW, LEVEL] = THRESHBONE(I, 'otsu') returns also the threshold LEVEL
%   calculated by the Otsu's method.
%
%   [BW] = THRESHBONE(I, 'adaptive', fsize, t, filterType, thresholdMode) 
%   Binarizes the input image I using the Wellner's Local Adaptive Threshold method.
%   fsize - Filter size used to determine the local weighted mean or local median.
%   t     - Depending on the value of 'mode' this is the value expressed as a percentage or fixed amount,
%               relative to the local average, or median  grey value, below which
%               the local threshold is set. Try values in the range -20 to +20. 
%   filterType     - Optional string specifying smoothing to be used: 'gaussian' or 'median'
%   thresholdMode  - Optional string specifying the way the threshold is
%               defined: 'relative' or 'fixed'
%
%   Class Support
%   -------------
%   BW is of class logical.
%
%   Example 1
%   -------
%      I = imread('cameraman.tif');
%      subplot(1,3,1);imshow(I);title('Original Image');
%      [BW, le] = threshbone(I, 'otsu');
%      subplot(1,3,2);imshow(BW);title('Otsu');
%      [BW, le] = threshbone(I, 'fixed', [140 180]);
%      subplot(1,3,3);imshow(BW);title('140 < GV < 180');
% 
%   Example 2
%   -------
%      I = imread('rice.png');
%      subplot(1,3,1);imshow(I);title('Original Image');
%      [BW, le] = threshbone(I, 'otsu');
%      subplot(1,3,2);imshow(BW);title('Otsu');
%      [BW] = threshbone(I, 'adaptive', 35, -35, 'median', 'fixed');
%      subplot(1,3,3);imshow(BW);title('LAT');
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   25/10/2016
%   Last update:  06/11/2018
%
%   adaptive method Copyright (c) 2008 Peter Kovesi
%   this function is part of the synchro toolbox    
%   ______________________________________________________

[method, p2, slicewise] = ParseInputs(varargin{:});

switch method
    case 'fixed'
        [BW, level] = fixedthresh(I, p2);
    case 'otsu'
        [BW, level] = otsuthresh(I, p2, slicewise);
    case 'ridler'
        [BW, level] = ridlerthresh(I, slicewise);
    case 'lakshmanan'
        [BW, level] = lakshmananthresh(I, p2);
    case 'adaptive'
        [BW] = kovesiadaptivethresh(I, p2);
end

function [BW, level] = otsuthresh(I, p2, slicewise)
    if slicewise
        % SLICEWISE Otsu's method
        fprintf('threshbone: SLICEWISE Otsu method\n');
        BW = false(size(I));
        for i=1:size(I,3)
            level = multithresh(I(:,:,i), p2);
            BW(:,:,i) = I(:,:,i) > level(end);
        end
    else
        % Otsu's method
        level = multithresh(I, p2);             % obtain threshold level with Otsu's Method
        fprintf('threshbone: Otsu method:\n');
        for i=1:length(level)
            fprintf('\tlevel %i: %2.2f\n', i, level(i));
        end
        BW = false(size(I));
        BW(I > level(end)) = true;
    end
end

function [BW, level] = lakshmananthresh(I, p2)
    % adaptive threshold for SAM images (Lakshmanan_2007)
    tentative_thresh = p2;                                               % [MRayl]
%     tentative_thresh = multithresh(I);                                               % [MRayl]
    I = double(I);
    Zmean_bone  = nanmean(I(I > tentative_thresh));                     % mean intensity of pixels exceeding tentative thresh
    Zstd_bone   = nanstd(I(I > tentative_thresh));                      % intensity standard deviation of pixels exceeding tentative thresh
    Zmean_pmma  = nanmean(I(I < tentative_thresh));                     % mean intensity of pixels below tentative thresh
    level = (Zmean_pmma + (Zmean_bone-Zstd_bone))/2;

    fprintf('threshbone: Adaptive method for SAM images (Lakshmanan_2007):\n\tlevel: %2.2f\n', level);
    BW = false(size(I));
    BW(I > level) = true;
end

function [BW, level] = ridlerthresh(I, slicewise)
    if slicewise
        % SLICEWISE ISODATA
        fprintf('threshbone: Iterative isodata SLICEWISE method\n');
        BW = false(size(I));
        for i=1:size(I,3)
            level = isodata(I(:,:,i));             % Compute global image threshold using iterative isodata method
            BW(:,:,i) = I(:,:,i) > level;
        end
    else
        % ISODATA
        level = isodata(I);             % Compute global image threshold using iterative isodata method
        fprintf('threshbone: Iterative isodata method:\n\tlevel: %2.2f\n', level);
        BW = false(size(I));
        BW(I > level) = true;
    end
end

function [BW, level] = fixedthresh(I, level)
    % Fixed global threshold (single or MAX+MIN level)
    if length(level) > 2     error('threshbone: Too Many Levels For This Method');

    elseif length(level) == 2
        % MIN and MAX global threshold
        fprintf('threshbone: 2-level global thresh:\n\tMIN: %2.2f\n\tMAX: %2.2f\n', min(level), max(level));
        mask1 = I > min(level);
        mask2 = I < max(level);
        BW = false(size(I));
        BW(mask1 & mask2) = true;

    elseif length(level) == 1
        % single valued global threshold
        fprintf('threshbone: single level global thresh:\n\tlevel: %2.2f\n', level);
        BW = false(size(I));
        BW(I > level) = true;
    end

end

function [BW] = kovesiadaptivethresh(I, p2)
    % Wellner's adaptive thresholding (Peter Kovesi)
    switch ndims(I)
        case 2
            fprintf('threshbone: Local Adaptive thresh:\n\t fsize: %i\n\t t: %i [%%]\n\t filterType: %s\n\t thresholdMode: %s\n', p2.fsize, p2.t, p2.filterType, p2.thresholdMode);
            [BW] = adaptivethresh(I, p2.fsize, p2.t, p2.filterType, p2.thresholdMode);
        case 3
            fprintf('threshbone: Slicewise Local Adaptive thresh:\n\t fsize: %i\n\t t: %i [%%]\n\t filterType: %s\n\t thresholdMode: %s\n', p2.fsize, p2.t, p2.filterType, p2.thresholdMode);
            BW = false(size(I));
            for i=1:size(I,3)
                BW(:,:,i) = adaptivethresh(I(:,:,i), p2.fsize, p2.t, p2.filterType, p2.thresholdMode);
            end
    end
end

function [method, p2, slicewise] = ParseInputs(varargin)

    % default values
    p2        = [];
    slicewise = false;

    % Check the number of input arguments.
    narginchk(0,6);

    % Determine threshold method from the user supplied string.
    if isempty(varargin)
        % fprintf('threshbone: No method given: I will use Otsu method.\n');
        method = 'otsu';
    else
        method = varargin{1};
        method = validatestring(method,{'fixed','otsu','lakshmanan','ridler','adaptive'},mfilename,'METHOD',1);
    end

    % default values
    switch method
        case 'otsu'
            p2 = 1;       % number of thresholds
        case 'lakshmanan'
            p2 = 4;       % [MRayl] tentative IMPEDANCE threshold
        case 'adaptive'
            p2.fsize = 3;
            p2.t = 15;
            p2.filterType = 'gaussian';
            p2.thresholdMode = 'relative';
    end

    switch nargin
        case 1
            % threshbone(I, 'otsu')
            % threshbone(I, 'lakshman')
            % threshbone(I, 'ridler')
            % threshbone(I, 'adaptive')

            switch method
                case {'fixed'}
                    error('threshbone: Not Enough Args For This Method')
            end

        case 2
            % threshbone(I, 'otsu', N)
            % threshbone(I, 'lakshman', tentative_thresh)
            % threshbone(I, 'fixed', level)
            % threshbone(I, 'adaptive', fsize)
            % threshbone(I, 'ridler', slicewise)
            
            p2 = varargin{2};
            
            switch method
                case {'lakshman'}
                    validateattributes(p2,{'numeric'},{'finite','real','nonempty'},...        %'scalar'
                        mfilename,'TENTATIVE_THRESH',3);
                case {'otsu'}
                    validateattributes(p2,{'double'},{'nonnegative','integer',...
                        'nonempty','finite'},...
                        mfilename,'N',3);
                    %   if  p2 > 1
                    %       error(message('images:fspecial:outOfRangeAlpha'))
                    %   end
                case {'fixed'}
                    validateattributes(p2,{'numeric'},{'finite','real','nonempty'},...        %'scalar'
                        mfilename,'LEVEL(S)',3);
                case {'adaptive'}
                    validateattributes(p2,{'char'},{'nonempty'},mfilename,'fsize',3);
                case {'ridler'}
                    slicewise = p2;
                    validateattributes(slicewise,{'logical'},{'nonempty'},mfilename,'slicewise',3);

            end
            
        case 3
            switch method
            % threshbone(I, 'otsu', N, slicewise)
            % threshbone(I, 'adaptive', fsize, t)
                case {'lakshman'}
                    error('threshbone: Too Many Args For This Method')
                case {'otsu'}
                    p2 = varargin{2};
                    slicewise = varargin{3};
                    validateattributes(p2,{'double'},{'nonnegative','integer',...
                        'nonempty','finite'},...
                        mfilename,'N',3);
                    validateattributes(slicewise,{'logical'},{'nonempty'},mfilename,'slicewise',3);
                case {'adaptive'}
                    p2.fsize = varargin{2};
                    p2.t = varargin{3};
                    validateattributes(p2.fsize,{'numeric'},{'nonempty'},mfilename,'fsize',3);
                    validateattributes(p2.t,{'numeric'},{'nonempty'},mfilename,'t',3);
            end
            
        case 4
            % threshbone(I, 'adaptive', fsize, t, filterType)
            p2.fsize = varargin{2};
            p2.t = varargin{3};
            p2.filterType = varargin{4};

        case 5
            % threshbone(I, 'adaptive', fsize, t, filterType, thresholdMode)
            p2.fsize = varargin{2};
            p2.t = varargin{3};
            p2.filterType = varargin{4};
            p2.thresholdMode = varargin{5};
            
    end

end

end
