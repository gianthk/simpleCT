function [Icrop, origin, cropsize] = crop(I, varargin)
%CROP (very basic) interactive crop
%   ICROP = CROP(I) interactively select crop region and crops input image I.
%
%   [ICROP, ORIGIN, CROPSIZE] = CROP(I) interactively select crop region and crops input image I.
%   Returns cropped data ICROP as well as ORIGIN and CROPSIZE.
%
%   [ICROP] = CROP(I, ORIGIN, CROPSIZE) Crops input image I from given ORIGIN and CROPSIZE.
%
%   Class Support
%   -------------
%   ICROP is of same class as I.
%
%   Example
%   -------
%      I = imread('cameraman.tif');
%      subplot(1,2,1);imshow(I);title('Original Image');
%      [Icrop] = crop(I);
%      subplot(1,2,2);imshow(Icrop);title('cropped');
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   29/10/2016
%   Last update:  09/03/2019
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

[origin, cropsize] = ParseInputs(varargin{:});

if isempty(origin)
    % run tool to select crop region
    [origin, cropsize] = selectcropregion(I);
end

switch ndims(I)
    case 2
        Icrop = I(origin(1):origin(1)+cropsize(1)-1, origin(2):origin(2)+cropsize(2)-1);
    case 3
        validateattributes(origin,{'double'},{'numel',3},mfilename,'ORIGIN',2);
        if origin(1)+cropsize(1)-1 > size(I,1),  error('crop: origin(1) + cropsize(1) - 1 > size(X,1)!');   end
        if origin(2)+cropsize(2)-1 > size(I,2),  error('crop: origin(2) + cropsize(2) - 1 > size(X,2)!');   end
        if origin(3)+cropsize(3)-1 > size(I,3),  error('crop: origin(3) + cropsize(3) - 1 > size(X,3)!');   end
        
        Icrop = I(origin(1):origin(1)+cropsize(1)-1, origin(2):origin(2)+cropsize(2)-1, origin(3):origin(3)+cropsize(3)-1);
end

function [origin, cropsize] = selectcropregion(I)
    switch ndims(I)
        case 2
            fprintf('select and double-click crop ROI!\n');
            figure;     imagesc(I); axis image
            title('select and double-click crop ROI')
            h = imrect;
            position = round(wait(h));
            y0 = position(1);   yd = position(3);
            x0 = position(2);   xd = position(4);
            z0 = [];            zd = [];
            close(gcf)
            clear h
            
            origin = [x0 y0];
            cropsize = [xn yn];

        case 3
            fprintf('select and double-click crop ROI!\n');
            figure;
%             subplot(1,3,1); imagesc(permute(max(I, [], 2), [3 1 2]));   axis image
%             subplot(1,3,2); imagesc(permute(max(I, [], 1), [3 2 1]));   axis image
            imagesc(max(I, [], 3)); axis image
            
            title('select and double-click crop ROI')
            h = imrect;
            position = round(wait(h));
            y0 = position(1);   yd = position(3);
            x0 = position(2);   xd = position(4);
            z0 = 1;             zd = size(I,3);
            close(gcf)
            clear h
            
            origin = [x0 y0 z0];
            cropsize = [xd yd zd];

    end


end

function [p1, p2] = ParseInputs(varargin)

    % default values
    p1        = [];
    p2        = [];

    % Check the number of input arguments.
    narginchk(0,2);
    
    switch nargin
        case 0
            % do nothing.. interactively select crop region

        case 1
            error('crop: Not Enough Args For This Method')

        case 2
            % origin and size are given
            p1 = varargin{1};
            p2 = varargin{2};
            
            validateattributes(p1,{'double'},{'nonnegative','nonzero','integer',...
                'nonempty','finite'},mfilename,'ORIGIN',2);

            validateattributes(p2,{'double'},{'nonnegative','nonzero','integer',...
                'nonempty','finite','numel',numel(p1)},mfilename,'CROPSIZE',3);
    end
    
end

end
