function [ CtPo ] = burghardt_CtPo( BW, C, voxelsize )
%BURGHARDT_CTPO calculates CtPo from bin + mask of cortical bone
%   the procedure is described in Burghardt et al. in 2010:
%         - Burghardt, Andrew J., et al. "Reproducibility of direct quantitative measures of cortical bone microarchitecture of the distal radius and tibia by HR-pQCT." Bone 47.3 (2010): 519-528.
%                 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2926164/
%
%   capital letters follow nomenclature of Fig. 1 in Burghardt et al.
% 
%   This procedure is also implemented in the Xcanco XCT IPL for cortical bone analysis.
%   Nevertheless, this function takes already binarized imgages (2D or 3D). 
%   You should keep in mind that the following steps are described in Burghardt et al.:
%       - filtering:    laplace-hamming in Burghardt et al.
%       - threshold:    40% of the max GV in Burghardt et al. (450 mgCA/cc in the Scanco IPL)
% 
%   voxelsize is assumed isotropic
%   ______________________________________________________
%
%   Author:         Gianluca Iori <gianthk.iori@gmail.com>
%   Created on:     07/04/2017
%   Last update:    26/04/2017
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________


%% ================ (3) CALCULATION OF CORTICAL MASK ================

D = BW & C;

% ================ (4) FIRST PORE ESTIMATE ================
%% 2D connectivity to select the cortical pores (removes pores connected to background)
a_thresh = 25;          % [mm]
a_thresh_pix = ceil(a_thresh/(voxelsize*voxelsize));

E = ~D;
for i = 1:size(D,3)
    Dinvslice = E(:,:,i);
    CC = bwconncomp(Dinvslice);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx = find(numPixels > a_thresh_pix);
    
    for j = 1:length(idx)
        Dinvslice(CC.PixelIdxList{idx(j)}) = 0;
    end
    E(:,:,i) = Dinvslice;
    
end
clear Dinvslice

%% ================ (5) SECOND PORE ESTIMATE ================
% cortical pores confidence map
poreconf = uint8(4*C - 3*D + 2*E);
poreconf(poreconf==2)=6;

%% hysteresis region growing in Z+
% starts from all pores assigned to first pore estimate and check for
% connectivity in Z+ to void voxels excluded from the first pore estimate
fprintf('slice:        res');
slices = size(C,3);
for k = 1:(slices-1)
    fprintf(1,[repmat('\b',1,numel(num2str(k))+numel(num2str(slices))) '\b%i/%i'], k, slices);
    poreconfslice = poreconf(:,:,k);
    [idr, idc] = find(poreconfslice == 6);
    for i = 1:length(idr)
        h = k;
        % search above each pixel included in first estimate. If a pixel
        % belongs to low confidence pore (value = 5) sets its value to 6
        row = idr(i);
        col = idc(i);
        while (h<slices && poreconf(row,col,h+1)==4)
            poreconf(row,col,h+1) = 5;
            % now that the low confidence pixel has been added to the first
            % pore estimate search in X-Y for connected neighbours with value=5
            poreconf(:,:,h+1) = poreconf(:,:,h+1) + uint8(grayconnected(poreconf(:,:,h+1),row,col,1));   % increase of 1 all pixels with intensity 4 < I < 6
            h = h+1;
        end
    end
end
fprintf('\n')
poreconf(poreconf>4)=6;

% now repeat the same procedure in Z-
fprintf('slice:        res');
slices = size(C,3);
for k = 1:(slices-1)
    fprintf(1,[repmat('\b',1,numel(num2str(k))+numel(num2str(slices))) '\b%i/%i'], k, slices);
    k0 = slices+1-k;
    poreconfslice = poreconf(:,:,k0);
    [idr, idc] = find(poreconfslice == 6);
    for i = 1:length(idr)
        h = k0;
        row = idr(i);
        col = idc(i);
        % search below each pixel included in first estimate. If a pixel
        % belongs to low confidence pore (value = 5) sets its value to 6
        while (h>1 && poreconf(row,col,h-1)==4)
            poreconf(row,col,h-1) = 5;
            % now that the low confidence pixel has been added to the first
            % pore estimate search in X-Y for connected neighbours with value=5
            poreconf(:,:,h-1) = poreconf(:,:,h-1) + uint8(grayconnected(poreconf(:,:,h-1),row,col,1));   % increase of 1 all pixels with intensity 6 < I < 8
            h = h-1;
        end
    end
end
fprintf('\n')
poreconf(poreconf>4)=6;

F = false(size(poreconf));
F(poreconf==6)=true;
clear poreconf poreconfslice

% ================ (6) FINAL POROSITY ESTIMATE ================
%% fill binary cortex with all found pores
G = F + D;

%% 2D connectivity to select small black pores (excludes external pores connected to background or marrow)
Ginv = ~G;
for i = 1:size(G,3)
    Ginvslice = Ginv(:,:,i);
    CC = bwconncomp(Ginvslice);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx = find(numPixels > a_thresh_pix);
    
    for j = 1:length(idx)
        Ginvslice(CC.PixelIdxList{idx(j)}) = 0;
    end
    Ginv(:,:,i) = Ginvslice;
    
end
clear Ginvslice G

%% add porosity from 2nd estimate
G = logical(Ginv + F);

%% discard pores with volume < 5 voxels
% H is to be used for Po.Dm and Po.Dm.SD calculation.. for simple Ct.Po
% estimates G is used
H = G;
RP = regionprops(H, 'Area');
CC = bwconncomp(H);
    
th = cat(1, RP.Area);
idx = find(th<5);

for j = 1:length(idx)
    H(CC.PixelIdxList{idx(j)}) = 0;
end

% ================ (7) CORTICAL MASK REFINEMENT ================
%% add pores to binary cortex for a refined cortical mask
% I will be used for quantitative measurements of Ct.Po and Ct.Th
I = logical(G + D);

%% ================ (8) Quantitative cortical bone analysis ================
% Intracortical pore volume (Ct.Po.V) is calculated as the integral volume 
% of all voxels identified as intracortical pore space (H), 
% while intracortical porosity (Ct.Po) is defined as a relative volumetric index of Ct.Po.V 
% normalized by the sum of the mineralized and intracortical pore volume.
TV = sum(sum(sum(I)));      % total volume
PV = sum(sum(sum(H)));      % pore volume

CtPo = (PV/TV)*100;

end

