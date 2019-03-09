function [dummy, bvtvout] = dummy(type, dummysize, voxelsize, poresize, bvtv, noisetype, noiseM, noiseS, blurS, graphics)
%DUMMY cerates in silico dummy of bone structure
%   [BW] = DUMMY('cortical') generates cortical bone structure dummy
%
%   [BW] = DUMMY('cortical', dummysize, voxelsize, poresize, bvtv, noisetype, noiseM, noiseS, blurS, graphics) 
%   dummysize -     output volume size
%   voxelsize -     (isotropic)
%   poresize -      diameter [mm]
%   bvtv -          target BV/TV [0-1]
%   noisetype
%   noiseM
%   noiseS
%
%   pippo = dummy('trabecular') generates trabecular bone (plate + rod-like) structure dummy
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   14/07/2017
%   Last update:  22/11/2017
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    noise = true;
    blur = true;
    if nargin < 10,         graphics = true;         	end
    if nargin < 9,
                            blur = false;
                            blurS = 0;
    end
    if isempty(blurS)
                            blur = false;
                            blurS = 0;
    end
    if nargin < 8,          noiseS = [];                end
    if nargin < 7,          noiseM = [];                end
    if nargin < 6,
                            noise = false;
                            noisetype = '';
                            noiseM = [];
                            noiseS = [];
    end
    if isempty(noisetype)   noise = false;              end
    if isempty(noiseM)      noiseM = 0;                 end
    if isempty(noiseS)      noiseS = 0.02;              end
    if nargin < 5,          bvtv = 0.9;                 end
    if isempty(bvtv),       bvtv = 0.9;                 end
    if nargin < 4,          poresize = 0.04;            end     % [mm]
    if isempty(poresize),   poresize = 0.04;            end
    if nargin < 3,          voxelsize = 0.01;        	end     % [mm]
    if isempty(voxelsize),  voxelsize = 0.01;           end
    if nargin < 2,          dummysize = [200 200 200];	end
    if isempty(dummysize),  dummysize = [200 200 200];  end
    if nargin < 1,          type = 'cortical';        	end
    if isempty(voxelsize)   voxelsize = 0.01;           end
    if isempty(dummysize)   dummysize = [200 200 200];  end
    if length(dummysize)==1 dummysize = [dummysize dummysize dummysize];    end
    
    switch type
        case 'cortical'
            % init dummy matrix
            dummy = true(dummysize);
            
            % create haversian canals until desired bvtv is reached
            while (sum(sum(sum(dummy)))/numel(dummy)) > bvtv
                % random Haversian Canal Center
                hcc(1) = ceil(rand(1)*dummysize(1));
                hcc(2) = ceil(rand(1)*dummysize(2));
                hr = poresize / voxelsize;        % healthy Haversian Canal: ~80 um
                % create pore (need KWave toolbox)
                tmp = dummy(:,:,1) & ~makeDisc(size(dummy,1), size(dummy,2), hcc(1), hcc(2), hr);
                dummy = repmat(tmp,1,1,size(dummy,3));
            end
            
            bvtvout = sum(sum(sum(dummy)))/numel(dummy);
            fprintf('BV/TV: %2.2f [%%]\n',100*bvtvout);
            
        case 'trabecular'
            dummy = false(dummysize);
            
            % create 1 rod + 1 plates
            hcc = [50 100];
            pc = [150];
            hr = 0.08 / voxelsize;        % Tb.Th: 160 um
            for i = 1:size(hcc,1)
                 for k = 1:size(dummy,3)
                     dummy(:,:,k) = dummy(:,:,k) | makeDisc(size(dummy,1), size(dummy,2), hcc(i,1), hcc(i,2), hr);
                     dummy(pc(i)-hr:pc(i)+hr,:,k) = true;
                 end
            end
            
            dummy = permute(dummy,[3,2,1]);
            
            % create 2 more rods + 2 plates
            hcc = [30 100; 80 100];
            pc = [130 170];
            hr = 0.08 / voxelsize;        % Tb.Th: 160 um
            for i = 1:length(hcc)
                 for k = 1:size(dummy,3)
                     dummy(:,:,k) = dummy(:,:,k) | makeDisc(size(dummy,1), size(dummy,2), hcc(i,1), hcc(i,2), hr);
                     dummy(pc(i)-hr:pc(i)+hr,:,k) = true;
                 end
            end
            
            bvtvout = sum(sum(sum(dummy)))/numel(dummy);
            fprintf('BV/TV: %2.2f [%%]\n',100*bvtvout);
    end
    
    %% blur
    if blur
        % gaussian blur with defined sigma
        dummy = imgaussfilt3(double(dummy), blurS);
    end
    
    %% add noise
    if noise
        dummy0 = double(dummy);
        dummy = imnoise(dummy0, noisetype, noiseM, noiseS);
        [PEAKSNR, SNR] = psnr(dummy, dummy0);
        fprintf('SNR = %2.2f \n',SNR);
    end
    
    %% plot
    if graphics
%         figure; slice(double(dummy),100,100,100)
        figure;
        subplot(1,3,1);
        imagesc(permute(dummy(:,round(size(dummy,2)/2),:),[1 3 2]))
        axis image
        colormap gray
        subplot(1,3,2);
        imagesc(permute(dummy(round(size(dummy,2)/2),:,:),[2 3 1]))
        axis image
        colormap gray
        subplot(1,3,3);
        imagesc(dummy(:,:,round(size(dummy,2)/2)))
        axis image
        colormap gray
    end
    
end
