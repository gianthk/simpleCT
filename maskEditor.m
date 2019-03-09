classdef maskEditor < handle
%MASKEDITOR manual edit of 3D masks
%   synchro class for manual editing of 3D masks:
%
%   The procedure is based upon displaying slices of a 3D dataset together
%   with the contour polygons of a given mask. The contour polygons are
%   manually editable.
%   The class works with an internal downsampled version of the mask to reduce
%   the number of slices that require manual correction. When the modified
%   mask is exported the downsampled version is interpolated to fit the
%   size of the original data.
%   Due to this, input data is cropped to a size which is multiple of the downsampling factor.
%   A GUI version of this procedure is also available.
%
%   This class is based on the use of the MATLAB's impoly.
%
%   See also MASKEDITOR_GUI, IMPOLY.
%   ______________________________________________________
%
%   Author: Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   03/05/2018
%   Last update:  09/05/2018
%
%   this class is part of the synchro toolbox    
%   ______________________________________________________

properties
                        % @type numeric
    im                  % input data
                        % 3D matrix of integers, doubles or logicals

                        % @type logical
    mask                % input mask
                        % 3D matrix of logicals; same size as the input data
                        
                        % @type logical
    isendo              % mask type: (0=peri; 1=endo)

                        % @type numeric
    fac                 % downsampling factor (default = 10)
    
                        % @type numeric
    rC                  % cropped rows

                        % @type numeric
    cC                  % cropped cols

                        % @type numeric
    sC                  % cropped slices

                        % @type numeric
    imC                 % cropped data
                        
                        % @type numeric
    sliceC              % slice for visualization
                        
                        % @type logical
    maskC               % cropped mask
                        
                        % @type logical
    maskCD              % cropped and downsampled mask

                        % @type logical
    maskCDe             % EDIT cropped and downsampled mask
    
                        % @type numeric
    contourCD           % contour coordinates of the cropped and downsampled mask
                        % [n x 3]
    
                        % @type numeric
    contourCDslice      % coordinates of 1 slice of contourCD
                        % [n x 2]
    
                        % @type handle
    poly                % polygone object for contour editing
                        
                        % @type numeric
    slicen              % the slice being edited
                        % integer

                        % @type numeric
    slicenD             % slice being edited in downsampled coordinates
                        % integer
    
                        % @type sctdata
    sctdataim           % sctdata class for handling data input and output

end

methods (Static = true, Access = private)
    
end

methods
    
    function cropFac(this)
        %CROPFAC mask -> maskC (crop data to fit downsampling factor)
        %   Crops input IM to a size which is multiple of the downsampling  factor FAC
%         if isempty(this.rC)
%             this.setRcsC;
%         end
        this.setRcsC;
        
        this.setImC( this.im(1:end-this.rC, 1:end-this.cC, 1:end-this.sC) );
        this.setMaskC( this.mask(1:end-this.rC, 1:end-this.cC, 1:end-this.sC) );
    end
    
    function padFac(this)
        %PADFAC maskC -> mask (adds back cropped voxels)
        mask = this.mask;
        mask(1:end-this.rC, 1:end-this.cC, 1:end-this.sC) = this.maskC;
        this.setMask(mask);
    end
    
    function calcAndSetContourCDslice( this )
        this.setContourCDslice(this.getContourCDslice(this.slicenD));
    end

    function initMaskCDe( this )
        this.maskCDe = [];
        for i=1:size(this.maskCD,3)
            this.maskCDe(:,:,i) = logical(imresize(single(this.maskCD(:,:,i)), size(this.maskC(:,:,1)))>0.5);
        end
    end

    function updateMaskCDe(this)
        if this.isendo
            this.maskCDe(:,:,this.slicenD) = ~logical(this.poly.createMask);
        else
            this.maskCDe(:,:,this.slicenD) = logical(this.poly.createMask);
        end
    end    

    function updateMaskCD(this)
        %UPDATEMASKCD maskCDe -> maskCD
        % slicewise resize to isotropic voxels
        if this.isendo
            for i=1:size(this.maskCDe,3)
                this.maskCD(:,:,i) = imresize(single(this.maskCDe(:,:,i)), 1/this.fac)>0.1;
            end
        else
            for i=1:size(this.maskCDe,3)
                this.maskCD(:,:,i) = imresize(single(this.maskCDe(:,:,i)), 1/this.fac)>0.5;
            end
        end
    end
    
    function updateMaskC(this)
        %UPDATEMASKC maskCD -> maskC
        % resize (back) the mask
        if this.isendo
            this.maskC = imresizen(single(this.maskCD), this.fac)>0.1;
        else
            this.maskC = imresizen(single(this.maskCD), this.fac)>0.5;
        end
    end

    function [fig, ca] = editContourCDslice(this)
        this.calcAndSetContourCDslice;
        this.setSliceC;
        figure;
        imagesc(this.sliceC);   colormap gray;  axis image;
        %imshowpair(this.sliceC, ~this.maskC(:,:,this.slicen));  axis image;
        fig = gcf;
        % set(cf,'units','normalized','outerposition',[0 0 1 1]);
        ca = gca;
        this.poly = impoly(ca, [this.fac*this.contourCDslice(:,1)-this.fac/2 this.fac*this.contourCDslice(:,2)-this.fac/2]);
    end

    function updateContourCDslice(this)
        this.contourCDslice = this.poly.getPosition;
    end
    
    function load(this)
        fprintf('Select data file..\n');
        this.importIm;
        currentdir = pwd;
        [path] = fileparts(this.sctdataim.headerfile);
        cd(path);
        fprintf('Select mask file..\n');
        this.importMask;
        cd(currentdir);
    end
    
    function importIm(this, filename)
        if nargin < 2
            this.sctdataim.load;
        else
            this.sctdataim.load(filename);
        end
        this.setIm(this.sctdataim.data);
    end
    
    function importMask(this, filename)
        if nargin < 2
            this.sctdataim.load;
        else
            this.sctdataim.load(filename);
        end
        this.setMask(logical(this.sctdataim.data));
    end
    
    function exportMask(this, filename)
        %EXPORTMASK maskCDe -> mask + EXPORT
        
        % from maskCDe to maskCD
        this.updateMaskCD;
        
        % from maksCD to maskC
        this.updateMaskC;
        
        % from maskC to mask
        this.padFac;
        
        % export
        this.sctdataim.setData(this.mask);
        if nargin < 2
            this.sctdataim.save;
        else
            this.sctdataim.save(filename);
        end
    end 

    % GETTERS   ___________________________________________________________
    function im = getIm(this)
        im = this.im;
    end
    
    function mask = getMask(this)
        mask = this.mask;
    end
    
    function fac = get.fac(this)
        fac = this.fac;
    end
    
    function rC = getRC(this)
        rC = this.rC;
    end
    
    function cC = getCC(this)
        cC = this.cC;
    end

    function sC = getSC(this)
        sC = this.sC;
    end
    
    function sliceC = getSliceC(this, sln)
    % returns volume slice at position (slice number) "sln" 
        validateattributes(sln,{'numeric'},{'nonempty','finite','positive'},mfilename,'sln',1);
        sliceC = this.imC(:,:,sln);
    end
    
    function [contourCD] = getContourCD(this)
        % gets the (cropped + downsampled) contour
        contourCD = this.contourCD;
    end
    
    function [contourCDslice] = getContourCDslice(this, slicenD)
        % gets 1 slice of the (cropped + downsampled) contour
        contour_idx = find(this.contourCD(:,3) == slicenD);
        contourCDslice = [this.contourCD(contour_idx, 2) this.contourCD(contour_idx, 1)];
    end
        
    % SETTERS   ___________________________________________________________
    function setIm( this, prop )
    % Setter for numeric property im
    % Usage:
    % 	 this.setIm(  prop )
    % Parameters:
    % 	this:  @type 
    % 	prop:  @type numeric
    
        validateattributes(prop,{'numeric', 'logical'},{'nonempty'},mfilename,'im');

        this.im = prop;
        fprintf('data set.\n');

%         this.setDims(ndims(prop));
%         datasize = size(prop);
%         this.setSize(datasize);
        
%         switch this.dims
%             case 3
%                 this.setXYslicen(round(datasize(3)/2));
%                 this.setXZslicen(round(datasize(1)/2));
%                 this.setYZslicen(round(datasize(2)/2));
%                 this.setOrthoslices;
% 
%             case 2
%                 this.setXYslicen(1);
%                 this.setXYslice;
%         end

    end
    
    function setMask( this, prop )
    % Setter for binary mask property
    % Usage:
    % 	 this.setMask(  prop )
    % Parameters:
    % 	this:  @type 
    % 	prop:  @type 
    
        validateattributes(prop,{'logical'},{'nonempty','finite'},mfilename,'mask');
        
        if size(prop) ~= size(this.im)
              error('maskEditor.setMask: mask must have same size as im!');
        end
        
        this.mask = prop;
        fprintf('mask set.\n');

        % check if the mask is endosteum or periosteum mask
        if prop(round(size(prop,1)/2), round(size(prop,2)/2), 1) == 1
            this.isendo = 0;
        else
            this.isendo = 1;
        end
        
        % create cropped versions to fit the downsampling factor
        this.cropFac;
        
        % after creating cropped version create also the downsampled one
        this.setMaskCD;

    end
    
    function setImC( this, prop )
        validateattributes(prop,{'numeric', 'logical'},{'nonempty'},mfilename,'imC');
        this.imC = prop;
        fprintf('cropped data set.\n');
        this.setSliceC;
    end
    
    function setMaskC( this, prop )
        validateattributes(prop,{'logical'},{'nonempty','finite'},mfilename,'maskC');
        this.maskC = prop;
        fprintf('cropped mask set.\n');
    end
    
    function setMaskCD( this )
        %SETMASKCD set cropped + downsampled mask
        % resize the mask and convert it again to logicals
        if this.isendo
            this.maskCD = imresizen(single(this.maskC), 1/this.fac) > 0.1;
        else
            this.maskCD = imresizen(single(this.maskC), 1/this.fac) > 0.5;
        end
        fprintf('cropped downsampled mask set.\n');
        % update the slice number in the downsampled coordinates
        this.setSlicenD(1);
        % draw the slice
        this.setSliceC;
        % initialize the editing mask
        this.initMaskCDe;
    end
    
    function setContourCD( this, prop )
    % Setter for contourCD property
    
        validateattributes(prop,{'numeric', 'logical'},{'nonempty'},mfilename,'contourCD');
        if size(prop,2) ~= ndims(this.im)
              error('maskEditor.setContourCD: contour must have 2 or 3 columns as the dimensions of the mask!');
        end
        this.contourCD = prop;
        fprintf('contour set.\n');
        this.calcAndSetContourCDslice;
    end
    
    function setContourCDslice( this, prop )
    % Setter for contourCDslice property
    
        validateattributes(prop,{'numeric', 'logical'},{},mfilename,'contourCDslice');     % ,{'nonempty'}
        if size(prop,2) ~= 2
              error('maskEditor.setContourCDslice: contour must have 2 columns!');
        end
        this.contourCDslice = prop;
        fprintf('contour slice set.\n');

    end
    
    function setEndosteumcontour ( this )
        % calculation of the mask (and contour) from the original size mask
        % fprintf('computing mask.. ');
        % tmp = endosteummask(this.maskC);
        % tmp = imresizen(single(tmp), 1/this.fac) > 0.5;
        % this.setContourCD(endosteumcontour(tmp));

        % periosteum contour calculation from the already downsampled mask
        this.setContourCD(endosteumcontour(this.maskCD));
        this.isendo = 1;
    end

    function setPeriosteumcontour ( this )
        % calculation of the mask (and contour) from the original size mask
        % fprintf('computing mask.. ');
        % tmp = periosteummask(this.maskC, 1, 1);
        % tmp = imresizen(single(tmp), 1/this.fac) > 0.5;
        % this.setContourCD(periosteumcontour(tmp));
        
        % periosteum contour calculation from the already downsampled mask
        this.setContourCD(periosteumcontour(this.maskCD));
        this.isendo = 0;
    end

    function setMaskCDe( this, prop )
        validateattributes(prop,{'logical'},{'nonempty','finite'},mfilename,'maskCDe');
        this.maskCDe = prop;
    end
    
    function setRC( this, prop )
        this.rC = prop;
    end

    function setCC( this, prop )
        this.cC = prop;
    end

    function setSC( this, prop )
        this.sC = prop;
    end

    function setRcsC( this )
        rcs = size(this.mask)-this.fac*floor(size(this.mask)/this.fac);
        this.setRC(rcs(1));
        this.setCC(rcs(2));
        this.setSC(rcs(3));
    end
    
    function setSlicen( this )
        this.slicen = round(this.slicenD*this.fac - this.fac/2);
    end
    
    function setSliceC( this )
        this.sliceC = this.getSliceC(this.slicen);
    end
    
    function setSlicenD( this, prop )
        this.slicenD = prop;
        validateattributes(prop,{'numeric'},{'nonempty', 'finite', 'positive'},mfilename,'slicenD');
        this.setSlicen;
    end
    
    function setFac( this, prop )
        this.fac = prop;
        this.cropFac;
        this.setMaskCD;
        this.initMaskCDe;
    end
    
    % CONSTRUCTOR   _______________________________________________________
    function this = maskEditor()
        this.im = [];
        this.mask = [];
        this.isendo = [];
        this.fac = 4;
        this.imC = [];
        this.sliceC = [];
        this.maskC = [];
        this.maskCD = [];
        this.maskCDe = [];
        this.contourCD = [];
        this.contourCDslice = [];
        this.poly = [];
        this.slicen = 10;
        this.slicenD = 1;
        switch exist('sctdata')
            case 2
                this.sctdataim = sctdata;
            otherwise
                warning('in order to use sctdata i/o methids the file sctdata.m should be in your MATLAB path.');
        end
    end
end

end