function listing = dir2(varargin)
%DIR2 list directory content without '.' and '..'
%   https://stackoverflow.com/questions/27337514/matlab-dir-without-and
%   Author:         Jubobs (https://stackoverflow.com/users/2541573/jubobs)
%   ______________________________________________________

    if nargin == 0
        name = '.';
    elseif nargin == 1
        name = varargin{1};
    else
        error('Too many input arguments.')
    end

    listing = dir(name);

    inds = [];
    n    = 0;
    k    = 1;

    while n < 2 && k <= length(listing)
        if any(strcmp(listing(k).name, {'.', '..'}))
            inds(end + 1) = k;
            n = n + 1;
        end
        k = k + 1;
    end

    listing(inds) = [];

end

