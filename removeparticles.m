function [BWout] = removeparticles(BW, varargin)
%REMOVEPARTICLES remove particles from binary image
%   BW2 = REMOVEPARTICLES(BW, METHOD) remove particles from input binary image BW with given METHOD.
%   Possible values for METHOD are:
%
%     'connectivity'  	(Default) remove particles based on the proportion of their size to their conenctivity to the image.
%                       Default value for connectivity is 10 (10 times more particle than contact pixels).
%
%     'size'            Remove particles based on their size. (Default size = 50 pixels). 
%
%   Depending on METHOD, REMOVEPARTICLES may take one additional parameter which you can supply.
%
%   BW2 = REMOVEPARTICLES(BW, 'connectivity', 20) remove particles with at least 20 times
%       more pixels than their pixels in contact with the rest of the image.
%
%   BW2 = REMOVEPARTICLES(BW, 'size', 10) remove particles with surface larger than 10 pixels.
%
%   Class Support
%   -------------
%   BW and BWout are of class logical.
%
%   Example 1
%   -------
%      BW = imread('cameraman.tif');
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
%   Created on:   08/11/2018
%   Last update:  08/11/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

[method, p2] = ParseInputs(varargin{:});

switch method
    case 'size'
        [BWout] = sizeparticles(BW, p2);
    case 'connectivity'
        [BWout] = connectivityparticles(BW, p2);
end

function [BWout] = sizeparticles(BW, p2)
    % remove H-connected pixels
    BW2 = bwmorph(BW, 'hbreak', 10);
    % remove small objects from binary image
    BW2 = bwareaopen(BW2,1000);
    BWout = BW2;
    % perform morphological opening (erosion followed by dilation) with object removal between
    BW2 = imerode(BW2, strel('square', 3));
    % remove small objects from binary image
    BW2 = bwareaopen(BW2,1000);
    % dilate back
    BW3 = imdilate(BW2, strel('square', 3));
    % tentative particles are the difference between BW3 nd BW
    particles = BW & not(BW3);
    % particle props
    props = regionprops(particles, 'Area', 'PixelIdxList');
    % delete particles with area exceeding given thresh
    for i=1:length(props)
        if props(i).Area > p2
            BWout(props(i).PixelIdxList) = false;
        end
    end
end

function [BWout] = connectivityparticles(BW, p2)
    % remove H-connected pixels
    BW2 = bwmorph(BW, 'hbreak', 10);
    % remove small objects from binary image
    BW2 = bwareaopen(BW2,1000);
    BWout = BW2;
    % perform morphological opening (erosion followed by dilation) with object removal between
    BW2 = imerode(BW2, strel('square', 2));
    % remove small objects from binary image
    BW2 = bwareaopen(BW2,1000);
    % dilate back
    BW3 = imdilate(BW2, strel('square', 2));
    % first guess particles are the difference between BW5 nd BW
    particles = BW & not(BW3);
    % particle props
    props = regionprops(particles, 'Area', 'PixelIdxList');
    % dilate the image without first guess particles
    BW3 = imdilate(BW3, strel('disk', 1));
    % contact pixels are pixels belonging to both particle and dilated BW6
    for i=1:length(props)
        ncontactpixels = sum(BW3(props(i).PixelIdxList));
        if props(i).Area/ncontactpixels > p2
            BWout(props(i).PixelIdxList) = false;
        end
    end
end

% _________________________________________________________________________
    
function [method, p2] = ParseInputs(varargin)

    % default values
    p2        = [];
    slicewise = false;

    % Check the number of input arguments.
    narginchk(0,2);

    % Determine threshold method from the user supplied string.
    if isempty(varargin)
        method = 'connectivity';
    else
        method = varargin{1};
        method = validatestring(method,{'size','connectivity'},mfilename,'METHOD',1);
    end

    % default values
    switch method
        case 'size'
            p2 = 50;        % [pixels] minimum particle size
        case 'connectivity'
            p2 = 10;        % 10 times more particle than contact pixels
    end

    switch nargin
        case 2
            p2 = varargin{2};
            
            switch method
                case {'size'}
                    validateattributes(p2,{'numeric'},{'finite','real','nonempty'},...        %'scalar'
                        mfilename,'PARTICLE_SIZE',3);
                case {'connectivity'}
                    validateattributes(p2,{'numeric'},{'finite','real','nonempty'},...        %'scalar'
                        mfilename,'PARTICLE_CONNECTIVITY',3);
            end
    end

end

end
