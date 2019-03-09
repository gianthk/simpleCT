function [ BWout, A ] = corticalmask( BW, voxelsize, erode3D, pad, dil_steps, th_thresh, dil_mm, a_thresh, A, verbose)
%CORTICALMASK calculates mask of cortical bone
%   [ BWout ] = CORTICALMASK( BW, voxelsize) mask of cortical compartment
%       from input binary image BW and (isotropic) vixelsize
%   [ BWout ] = CORTICALMASK( BW, voxelsize, 1) activate 3D erosion
%       with pseudosphere structuring element (default = 2D; 'disk')
%   [ BWout ] = CORTICALMASK( BW, voxelsize, [], 1) activate image zeros-padding (default = off)
%   [ BWout ] = CORTICALMASK( BW, voxelsize, [], [], 3) control number of dilate/erosion sub-steps (default = 5) 
%   [ BWout ] = CORTICALMASK( BW, voxelsize, [], [], [], 0.22) control thickness threshold for filling thin pores (default = 0.164 mm [Burghardt et. al. 2010])
%   [ BWout ] = CORTICALMASK( BW, voxelsize, [], [], [], [], 1.5) control thickness of dilate/erode element (default = 2.46 mm [Burghardt et. al. 2010])
%   [ BWout ] = CORTICALMASK( BW, voxelsize, [], [], [], [], [], 10) control minimum square area of unconnected structs to be removed (default = 25 mm^2)
%   [ BWout ] = CORTICALMASK( BW, voxelsize, [], [], [], [], [], [], A) bypass calculation of periosteum mask (A) with the Burghardt procedure and feed it from the user
%   [ BWout ] = CORTICALMASK( BW, voxelsize, [], [], [], [], [], [], [], 1) activate graphical output (default = off)
%   [ BWout, A ] = CORTICALMASK( BW, voxelsize, ... ) output periosteal mask
%
%   the procedure is described in Burghardt et. al. in 2010:
%         - Burghardt, Andrew J., et al. "Reproducibility of direct quantitative measures of cortical bone microarchitecture of the distal radius and tibia by HR-pQCT." Bone 47.3 (2010): 519-528.
%                 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2926164/
%
%   capital letters follow nomenclature of Fig. 1 in Burghardt et. al.
% 
%   This procedure is also implemented in the Xcanco XCT cortical analysis routine.
%   Nevertheless, this function takes already binarized imgages (2D or 3D). 
%   You should keep in mind that the following steps are prescribed by Burghardt et. al.:
%       - filtering:    laplace-hamming in Burghardt et. al.
%       - threshold:    40% of the max GV in Burghardt et. al., 450 mgCA/cc in the Scanco routine
% 
%   voxelsize is assumed isotropic
%   ______________________________________________________
%
%   Author:         Gianluca Iori <gianthk.iori@gmail.com>
%   Created on:     07/04/2017
%   Last update:    09/03/2019
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    if nargin<10                verbose = 0;                        end
    if nargin<9                 A = [];                             end
    if nargin<8                 a_thresh = 25;                      end     % [mm^2] area threshold for removing thin sparse struts
    if isempty(a_thresh)        a_thresh = 25;                      end
    if nargin<7                 dil_mm = [2.46, 2.46];              end     % [mm] imdilate 2.46 mm (Burghardt et. al. 2010)
    if isempty(dil_mm)          dil_mm = [2.46, 2.46];              end
    if numel(dil_mm) == 1,      dil_mm = [dil_mm, dil_mm];          end
    if nargin<6                 th_thresh = 0.164;                  end
    if isempty(th_thresh)       th_thresh = 0.164;                  end
    if nargin<5                 dil_steps = [10 10];                end
    if isempty(dil_steps)       dil_steps = [10 10];                end
    if numel(dil_steps) == 1,   dil_steps = [dil_steps, dil_steps]; end
    if nargin<4                 pad = false;                        end
    if isempty(pad)             pad = false;                        end
    if nargin<3                 erode3D = false;                    end
    if isempty(erode3D)         erode3D = false;                    end

    % constants
    dil_voxels = round((dil_mm/voxelsize)/2);

    % add padding
    if pad
        padsize = max(dil_voxels) + 10;
        BW = padarray(BW,[padsize padsize 0]);
        if verbose>0    fprintf('added padding \n');    end
    end

    %% ================ (0) CLEANING ================
    % remove sparse structs around mask
    if verbose>0
        fprintf('(0) CLEANING \nremove sparse structs around mask... \n');
    end
    a_thresh_pix = ceil(a_thresh/(voxelsize*voxelsize));

    for i = 1:size(BW,3)
        BWslice = BW(:,:,i);
        CC = bwconncomp(BWslice);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        idx = find(numPixels < a_thresh_pix);

        for j = 1:length(idx)
            BWslice(CC.PixelIdxList{idx(j)}) = 0;
        end
        BW(:,:,i) = BWslice;

    end
    clear BWslice
    
    if isempty(A)
        %% ================ (1) EXTRACTION OF PERIOSTEAL SURFACE ================
        if verbose>0    fprintf('(1) EXTRACTION OF PERIOSTEAL SURFACE \n');     end

        % dilate bone to remove Volksmann's canals
        
        if verbose>0
            fprintf('dilate bone to remove Volksmann''s canals... ');
            tic()
        end

        if (erode3D && ndims(BW)==3)
            se_voxels = round(dil_voxels(1)/dil_steps(1));
            se = strel3d(2*se_voxels);           % 3D imdilate (EXTREMELY slow)
            A = BW;
            if verbose>0
                fprintf('(3D)  ');
            end

            for i=1:dil_steps(1)
                A = imdilate(A, se);
            end
        else
            se = strel('disk', dil_voxels(1));
            A = imdilate(BW, se);
        end    
        if verbose>0
            toc()
        end   

        % 2D connectivity to select the background (fills the bone; keeps only background)
        Ainv = ~A;
        A = true(size(Ainv));
        for i = 1:size(A,3)
            CC = bwconncomp(Ainv(:,:,i));
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [biggest,idx] = max(numPixels);
            Aslice = true(size(Ainv,1),size(Ainv,2));
            Aslice(CC.PixelIdxList{idx}) = 0;
            A(:,:,i) = Aslice;
        end
        clear Ainv Aslice

        % erode back 2.46 mm to original size
        if verbose>0
            fprintf('erode back to original size... ');
            tic()
        end
        if (erode3D && ndims(BW)==3)
            if verbose>0    fprintf('(3D)  ');  end
            for i=1:dil_steps(1)
                A = imerode(A, se);
            end
        else
            A = imerode(A, se);
        end
        if verbose>0
            toc()
        end
    else
        if verbose>0    fprintf('..extraction of periosteal surface (1) skipped.. \n');  end
        % add padding
        if pad
            A = padarray(A,[100 100 0]);
            if verbose>0    fprintf('added padding to periosteal surface\n');    end
        end
    end

    %% ================ (2) EXTRACTION OF ENDOSTEAL SURFACE ================
    if verbose>0    fprintf('(2) EXTRACTION OF ENDOSTEAL SURFACE \n');  end

    % mask binary image with A
    B = ~BW.*A;

    % 2D connectivity to select the background (removes trabeculae)
    if verbose>0    fprintf('2D connectivity to select the background (removes trabeculae) \n');    end
    Binv = ~B;
    B = true(size(Binv));
    for i = 1:size(B,3)
        CC = bwconncomp(Binv(:,:,i));
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        Bslice = true(size(Binv,1),size(Binv,2));
        Bslice(CC.PixelIdxList{idx}) = 0;
        B(:,:,i) = Bslice;
    end

    clear Binv Bslice

    % calculate thickness and subtract thin elements
    if verbose>0    fprintf('calculate thickness and subtract thin elements \n');    end
    th_thresh_pix = th_thresh / voxelsize;

    for i = 1:size(B,3)
        Bslice = B(:,:,i);
        CC = bwconncomp(Bslice);
        RP = regionprops(Bslice, 'MinorAxisLength');
        th = cat(1, RP.MinorAxisLength);
        idx = find(th<th_thresh_pix);

        for j = 1:length(idx)
            Bslice(CC.PixelIdxList{idx(j)}) = 0;
        end
        B(:,:,i) = Bslice;
    end
    clear Bslice

    % remove unconnected pores
    B = removenonconnectedstruts(imopen(B, strel3d(5)), 'regionprops');
    
    % dilate marrow to connect
    if verbose>0    fprintf('dilate marrow to connect \n');    end
    if (erode3D && ndims(BW)==3)
        se_voxels = round(dil_voxels(2)/dil_steps(2));
        se = strel3d(2*se_voxels);           % 3D imdilate (slow)
        if verbose>0
            fprintf('dilate 3D endo... ');
            tic()
        end

        for i=1:dil_steps(2)
            B = imdilate(B, se);
        end
        if verbose>0
            toc()
        end
    else
        se = strel('disk', dil_voxels(2));
        if verbose>0
            fprintf('dilate 2D endo... ');
            tic()
        end    
        B = imdilate(B, se);
        if verbose>0
            toc()
        end
    end

    % 2D connectivity to select the background (removes holes)
    if verbose>0    fprintf('2D connectivity to select the background (removes holes) \n');    end
    Binv = ~B;
    B = true(size(Binv));

    for i = 1:size(B,3)
        CC = bwconncomp(Binv(:,:,i));
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [biggest,idx] = max(numPixels);
        Bslice = true(size(Binv,1),size(Binv,2));
        Bslice(CC.PixelIdxList{idx}) = 0;
        B(:,:,i) = Bslice;
    end

    clear Binv Bslice

    % erode back to original size
    if erode3D
        if verbose>0
            fprintf('erode back to original size (3D)...');
            tic()
        end
        for i=1:dil_steps(2)
            B = imerode(B, se);
        end
        if verbose>0
            toc()
        end
    else
        if verbose>0    fprintf('erode back to original size (2D)...');
            tic()
        end
        B = imerode(B, se);
        if verbose>0
            toc()
        end
    end

    %% ================ (3) CALCULATION OF CORTICAL MASK ================
    if verbose>0    fprintf('(3) CALCULATION OF CORTICAL MASK.. ');        end
    % PERI - ENDO
    BWout = A & ~B;

    % remove padding
    if pad
        BWout = BWout(padsize+1:end-padsize,padsize+1:end-padsize,:);
        A = A(padsize+1:end-padsize,padsize+1:end-padsize,:);
    end

    if verbose>0    fprintf('Done!\n');        end

    %% ================ (4) FINAL CLEANING ================
    % remove sparse structs around mask
    if verbose>0
        fprintf('(4) FINAL CLEANING \nremove sparse structs around mask... \n');
    end
    
    for i = 1:size(BWout,3)
        BWslice = BWout(:,:,i);
        CC = bwconncomp(BWslice);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        idx = find(numPixels < a_thresh_pix);

        for j = 1:length(idx)
            BWslice(CC.PixelIdxList{idx(j)}) = 0;
        end
        BWout(:,:,i) = BWslice;
    end
    
    %% visualization video
    %     figure;
    %     for i = 1:size(BW,3)
    %     %     imagesc(I(:,:,i)+C(:,:,i)+5*G(:,:,i));
    % %         imagesc(BW22(:,:,i)+pippo.mask(:,:,i));
    %         imagesc(A(:,:,i));
    %         pause(0.1);
    %     end

    %     for sln=[70:110] figure; imagesc(BW(:,:,sln)+BWout(:,:,sln)); colormap gray; title('BW'); end

end

