function [ctth, ctth_array, ctth_bins, ctth_ksdensity, method] = CtTh(cloudA, cloudB, varargin)
%CTTH Cortical Thickness
%   ctth = CtTh(cloudA, cloudB, METHOD) calculates cortical bone thickness (Ct.Th)
%   between points cloudA and cloudB (endosteum and periosteum) using specified METHOD.
%   cloudA and cloudB are defined as nx2 or nx3 coordinates arrays.
%   Possible values for METHOD are:
%
%     'pdist2'      (Default) Euler distance between each cloudA-cloudB pair of points.
%                   The minimum distance of cloudB to each cloudA point is
%                   stored. CtTh is the maximum occurrance of minimum dist [1].
%     'cloudnorm'   normals to the points cloud.
%     'fitsphere'   
%
%   Depending on METHOD, CTTH may take one additional parameter
%   which you can supply.
%
%   [ctth] = CtTh(cloudA, cloudB, METHOD, 1) activates graphical output.
%
%   [ctth] = CtTh(cloudA, cloudB, 'cloudnorm', 0, nsamples, tol) calculates thickness
%   by projecting a normal from nsamples points of cloudA onto cloudB. The
%   cloudA normals are calculated by plane-fitting the point cloud within
%   tol distance.
%
%   Class Support
%   -------------
%   cloudA and cloudB are converted to doubles and should have 2 or 3 columns (2D or 3D).
%
%   Example 1: mouse femur 3D Ct.Th
%   -------
%      data = MHDdata.load('mouse_femur2.mhd');
%      data = data(1:40,1:60,60:100);
%      mask = threshbone(data);
%      peri = periosteumcontour(mask);
%      endo = endosteumcontour(mask);
%      
%      bdat = [20 4 10; 20 35 10; 35 35 10; 35 4 30; 20 4 30; 20 35 30; 35 35 30; 35 4 30];
%      in_endo = inpolygon(endo(:,1), endo(:,2), bdat(:,2), bdat(:,1));
%      in_peri = inpolygon(peri(:,1), peri(:,2), bdat(:,2), bdat(:,1));
%      
%      CtTh(endo(in_endo,:), peri(in_peri,:), 'pdist2', 1)
%
%   Example 2: 2D Ct.Th distribution from SAM image
%   -------
%      mask = imread('S:\AG\AG-Raum\Daten\Tacosound\SAM\final_aligned\femur_shaft\masks\cort\1955_L.tif');
%      peri = periosteumcontour(mask);
%      endo = endosteumcontour(mask);
%      CtTh(endo, peri, 'pdist2', 1)
%
%   References
%   -------
%      1.   Chappard, C., S. Bensalah, C. Olivier, P. J. Gouttenoire, A. Marchadier, C. Benhamou, and F. Peyrin.
%           “3D Characterization of Pores in the Cortical Bone of Human Femur in the Elderly at Different Locations as Determined by Synchrotron Micro-Computed Tomography Images.”
%           Osteoporosis International 24, no. 3 (March 1, 2013): 1023–33.
%           https://doi.org/10.1007/s00198-012-2044-4.
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   25/04/2018
%   Last update:  19/12/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

if isempty(cloudA) || isempty(cloudB)
    ctth = NaN;
    ctth_array = NaN;
    return
end

[method, p2] = ParseInputs(varargin{:});

switch method
    case 'pdist2'
        [ctth, ctth_array, ctth_bins, ctth_ksdensity] = pdist2th(cloudA, cloudB, p2.verbose, p2.parameter);
    case 'cloudnorm'
        warning('CtTh: cloudnorm method not implemented yet..');
        return
        [ctth, ctth_array] = cloudnormth(cloudA, cloudB, p2);
    case 'fitsphere'
        warning('CtTh: fitsphere method not implemented yet..');
        return
        [ctth, ctth_array] = fitsphere(cloudA, cloudB, p2);
end

function [ctth, ctth_array, ctth_bins, ctth_ksdensity] = pdist2th(cloudA, cloudB, verbose, parameter)
    if nargin < 4,          parameter = 'peak';    end
    if isempty(parameter),  parameter = 'peak';    end
    if nargin < 3,          verbose = 0;    end
    if isempty(verbose),    verbose = 0;    end
    
    % distance method
    % calculate distance matrix
    dist_matrix = pdist2(cloudA, cloudB);

    % mimimum distances from each cloudA point
    ctth_array = min(dist_matrix');
    
    % CtTh ad the mimimum distance with most occurrances (ksdensity peak)
    [ctth_ksdensity, ctth_bins, bw] = ksdensity(ctth_array);
    npeaks = numel(findpeaks(ctth_ksdensity));
    ctth_peaks = sort(findpeaks(ctth_ksdensity), 'descend');
    
    % if more than one peak exist and the most common peak isn't at least twice the second 
    while npeaks>1 && ctth_peaks(2)>0.5*ctth_peaks(1) 
        % 10% BandWidth increase
        bw = 1.1*bw;
        [ctth_ksdensity, ctth_bins, bw] = ksdensity(ctth_array, 'Bandwidth', bw);
        npeaks = numel(local_max(ctth_ksdensity));
    end
    
    % [ctth_hist, ctth_bins] = histcount(ctth_array);
    peak_idx = find(ctth_ksdensity == max(ctth_ksdensity));
    % peak_idx = find(min(abs(ctth_ksdensity-max(ctth_ksdensity))) == abs(ctth_ksdensity-max(ctth_ksdensity)));
    
    switch parameter
        case 'peak'
            ctth = ctth_bins(peak_idx);
        case 'median'
            ctth = median(ctth_array);
        case 'mean'
            ctth = mean(ctth_array);
        case 'min'
            ctth = min(ctth_array);
    end
    
    % graphical output (optional)
    if verbose > 0
        figure;
        subplot 121;
        switch size(cloudA, 2)
            case 2
                plot(cloudA(:,1), cloudA(:,2), '.b', cloudB(:,1), cloudB(:,2), '.r');
                axis image;
                % set(gca,'Ydir','Reverse')
            case 3
                plot3(cloudA(:,1), cloudA(:,2), cloudA(:,3), '.b')
                hold on;
                plot3(cloudB(:,1), cloudB(:,2), cloudB(:,3), '.r')
        end
        axis image;
        subplot 122;
        plot(ctth_bins, ctth_ksdensity, '.-k', 'LineWidth', 1.5);
        hold on;
        plot(ctth, max(ctth_ksdensity), 'or', 'LineWidth', 1.5);
        plot([ctth ctth], [0 max(ctth_ksdensity)], 'r', 'LineWidth', 1.5);
        xlabel('minimum thickness [pixels]');
        ylabel('frequency');        
    end
    
end
    
function [ctth, ctth_array] = cloudnormth(cloudA, cloudB, verbose, nsamples, tol);
    
end

% _________________________________________________________________________
    
function [method, p2] = ParseInputs(varargin)

    % default values
    p2        = [];

    % Check the number of input arguments.
    narginchk(0,4);

    % Determine method from the user supplied string.
    if isempty(varargin)
        method = 'pdist2';
    else
        method = varargin{1};
        if isempty(method), method = 'pdist2';  end
        method = validatestring(method,{'pdist2','cloudnorm','spherefit'},mfilename,'METHOD',1);
    end

    % default values
    switch method
        case 'pdist2'
            p2.verbose = 0;
            p2.parameter = 'peak';
        
        case 'cloudnorm'
            p2.verbose = 0;
            p2.nsamples = [];
            p2.tol = [];
    end

    switch nargin
        case 1
            % CtTh(I, 'pdist2')
            % CtTh(I, 'cloudnorm')
            % CtTh(I, 'spherefit')

        case 2
            % ACTIVATE GRAPHICAL OUTPUT
            % CtTh(I, 'pdist2', 1)
            % CtTh(I, 'cloudnorm', 1)
            % CtTh(I, 'spherefit', 1)

           p2.verbose = varargin{2};
           p2.parameter = '';

           switch method
              case {'pdist2'}
                  validateattributes(p2.verbose,{'numeric'},{'nonnegative','integer','nonempty','finite'},mfilename,'VERBOSE',3);
              case {'cloudnorm'}
                  validateattributes(p2.verbose,{'numeric'},{'nonnegative','integer','nonempty','finite'},mfilename,'VERBOSE',3);
              case {'spherefit'}
                  validateattributes(p2.verbose,{'numeric'},{'nonnegative','integer','nonempty','finite'},mfilename,'VERBOSE',3);
           end
           
        case 3
            % CtTh(I, 'pdist2', [], 'median')
            % CtTh(I, 'cloudnorm', [], nsamples)
            % CtTh(I, 'spherefit', [], )

           switch method
              case {'pdist2'}
                  % error('CtTh: Too Many Args For This Method')
                  p2.verbose = varargin{2};
                  p2.parameter = varargin{3};
                  
                  validateattributes(p2.verbose,{'double'},{'nonnegative','integer','nonempty','finite'},mfilename,'VERBOSE',3);
                  validateattributes(p2.parameter,{'char'},{'nonempty'},mfilename,'PARAMETER',4);
                  
                  p2.parameter = validatestring(p2.parameter,{'peak','median','mean','min'},mfilename,'PARAMETER',1);
              case {'cloudnorm'}
                  p2.verbose = varargin{2};
                  p2.nsamples = varargin{3};
                  
                  validateattributes(p2.verbose,{'double'},{'nonnegative','integer','nonempty','finite'},mfilename,'VERBOSE',3);
                  validateattributes(p2.nsamples,{'double'},{'nonnegative','integer','nonempty','finite'},mfilename,'NSAMPLES',4);
              case {'spherefit'}
                  % validateattributes(p2,{'numeric'},{'finite','real','nonempty'},mfilename,'LEVEL(S)',3);
           end
           
        case 4
            % CtTh(I, 'cloudnorm', [], nsamples, tol)
            p2.verbose = varargin{2};
            p2.nsamples = varargin{3};
            p2.tol = varargin{4};
            switch method
              case {'pdist2'}
                  error('CtTh: Too Many Args For This Method')
              case {'spherefit'}
                  error('CtTh: Too Many Args For This Method')
            end

    end

end

end
