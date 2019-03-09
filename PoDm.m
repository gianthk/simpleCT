function [dm_array] = PoDm(BW, verbose)
%PoDm Pore Diameter
%   podm = PoDm(BW)     returns array of minimum pore diameters given a binary image of the pores
%   podm = PoDm(BW, 1)  activate graphical output
%
%   Class Support
%   -------------
%   BW should be a 2D matrix of logicals.
%
%   Example 1: trabecular thickness (beta)
%   -------
%       r = max_inscribed_circle(bwperim(~threshbone(imread('spongy.png'))),1);
%
%   Example 2: Po.Dm array; comparison with Equivalent Diamater from Matlab's regioprops
%   -------
%       pores = imread('fivepores.tif');
%       stats = regionprops(pores, 'EquivDiameter');
%       for i=1:length(stats), PoDm_equivdiam(i)=stats(i).EquivDiameter; end
%       PoDm_circle = PoDm(pores, 1);
%       regressionplot(PoDm_equivdiam, PoDm_circle);
%
%   References
%   -------
%      1.   Maximum Inscribed Circle using Distance Transform   version 1.0.0.0 (12.6 KB) by Tolga Birdal
%           https://de.mathworks.com/matlabcentral/fileexchange/30805-maximum-inscribed-circle-using-distance-transform
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   25/07/2018
%   Last update:  25/07/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

if nargin<2,            verbose = false;    end
if isempty(verbose),    verbose = false;    end

PoSTATS = regionprops(BW, 'MinorAxisLength', 'Image');
% figure; imagesc(PoSTATS(103).Image)

dm_array = zeros(length(PoSTATS),1);

if verbose
    % graphical output ON
    for i=1:length(PoSTATS)
    %     BW2 = PoSTATS(i).Image;
        BW2 = imresize(PoSTATS(i).Image,4);
        BW2p = bwperim(BW2);
        if all(BW2(:) == BW2p(:))
            dm_array(i) = PoSTATS(i).MinorAxisLength;
        else
    %         dm_array(i) = 2*max_inscribed_circle(BW2p);
            dm_array(i) = 0.5*max_inscribed_circle(BW2p, 1);
        end
    end
else
    % graphical output OFF
    for i=1:length(PoSTATS)
    %     BW2 = PoSTATS(i).Image;
        BW2 = imresize(PoSTATS(i).Image,4);
        BW2p = bwperim(BW2);
        if all(BW2(:) == BW2p(:))
            dm_array(i) = PoSTATS(i).MinorAxisLength;
        else
    %         dm_array(i) = 2*max_inscribed_circle(BW2p);
            dm_array(i) = 0.5*max_inscribed_circle(BW2p, 0);
        end
    end
end

dm_array = double(dm_array);

end
