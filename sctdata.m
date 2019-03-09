classdef sctdata < handle
%   sctdata class for 2D and 3D computed tomography (CT) data
%
%   example call:   tmp = sctdata       % initialize sctdata object
%   ______________________________________________________
%
%   Author:         Gianluca Iori <gianthk.iori@gmail.com>
%   BSRT - Charite Berlin
%   Created on:   13/10/2016
%   Last update:  09/03/2019
%
%   this class is part of the synchro toolbox    
%   ______________________________________________________

properties
    
    % @type numeric
    data                % 2D or 3D matrix of pixeldata
                        % [rows columns slices] of the matrix correspond to
                        % [Y X Z] dimensions of the dataset
    
    % @type logical
    mask                % 2D or 3D matrix of binary data
                        % [rows columns slices] of the matrix correspond to
                        % [Y X Z] dimensions of the dataset

    % @type numeric
    min_data_value      % minimum pixel value
    
    % @type numeric
    max_data_value      % minimum pixel value

    % @type numeric
    XYslice             % 2D XY-slice through volume data
    
    % @type numeric
    XZslice             % 2D XZ-slice through volume data
    
    % @type numeric
    YZslice             % 2D XZ-slice through volume data
    
    % @type numeric
    XYslicemask         % 2D XY mask slice through volume data
    
    % @type numeric
    XZslicemask         % 2D XZ mask slice through volume data
    
    % @type numeric
    YZslicemask         % 2D XZ mask slice through volume data
    
    % @type numeric
    XYslicen            % XY-slice number
    
    % @type numeric
    XZslicen            % XZ-slice number
    
    % @type numeric
    YZslicen            % YZ-slice number
    
    % @type numeric
    min_xproj           % 2D image of minimum X-projection of dataset pixel data

    % @type numeric
    max_xproj           % 2D image of maximum X-projection of dataset pixel data
    
    % @type numeric
    min_yproj           % 2D image of minimum Y-projection of dataset pixel data

    % @type numeric
    max_yproj           % 2D image of maximum Y-projection of dataset pixel data
    
    % @type numeric
    min_zproj           % 2D image of minimum Z-projection of dataset pixel data

    % @type numeric
    max_zproj           % 2D image of maximum Z-projection of dataset pixel data
    
    % @type numeric
    dims                % scalar representing number of data dimensions
    
    % @type numeric
    size                % integer array of [X Y Z] data size ([ncolumns nrows nslices] of data)
    
    % @type numeric
    voxelsize           % double array of [X Y Z] voxel size
    
    % @type char
    samplename = '';    % sample name
    
    % @type char
    elementtype = '';   % data type

    % @type char
    headerfile = '';    % source file location and filename
    
    % @type char
    rawfile = '';       % source file (RAW, if different than headerfile) location and filename
    
    % @type numeric
    transformmatrix;    % 9x1 array containing 3D data rotation matrix
                        % for an homogeneous coordinates system
                        % tipically: TransformMatrix = [1 0 0 0 1 0 0 0 1] (NO ROTATION)
                        % see MHD header format for more info
    
    % @type numeric
    centerofrotation;   % 3x1 array of center of rotation coordinates
                        % tipically: centerofrotation = [0 0 0]
    
    % @type numeric
    offset              % 3x1 array of center of coordinates of starting data pixels
                        % tipically: offset = [0 0 0]
    
    % @type numeric
    dataoffset          % file bytes offset to start of data
    
    % @type numeric
    intensity           % [uA]
    
    % @type numeric
    energy              % [kV]
    
    % @type numeric
    mu_scaling          % -------------------------------------------------
                        % these 3 parameters allow conversion of Hounsfield Units
    
    % @type numeric
    slope_dens          % into BMD data. the calibration rule (valid for scanco
                        % scanners and if a phantom calibration was performed) is:
    
    % @type numeric
    offset_dens         % BMD = HU*(slope_dens/muScaling) + offset_dens;
                        
    % @type numeric
    nr_of_projections   % int
    
    % @type numeric
    nr_of_bytes         % int
    
    % @type logical
    binarydata          % (true if data is binary)
    
    % @type numeric
    thresholdvalue      % records value used for thresholding (if mask is present)
    
    % @type logical
    BMDdata             % (true if data calibrated to BMD)
    
end

properties (Access = private)
    pathstr
    name
    ext
end

properties (Dependent)
%         C
end

methods
    %% sctdata modify methods
    function [cropobj] = crop(this, x0, y0, z0, xd, yd, zd)
    %cropobj = crop(this, x0, y0, z0, xd, yd, zd)
    %   sct crop method
    %
    %   WARNING: offset not implemented yet.. origin information will be lost
    %
    %   inputs are the first crop point (x0, y0, z0) and size of the crop region (xd, yd, zd), in mm
    %       if output argument is given crop returns sctdata
    %       object of the cropped portion. If no output argument is given
    %       crop overwrites the input sctdata object.
    %
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   pippo.crop(40, 60, 50, 30, 30, 20).showslice
    %                       crops a volume of 30x30x20 [mm] from the point of
    %                       coordinates (40, 60, 50) [mm]
    %
    %                   pippo.crop      launches mini GUI for interactive
    %                                   selection of crop region
    %   ______________________________________________________
    %
    %   Created on:   29/10/2016
    %   Last update:  11/01/2018
    %   ______________________________________________________
    
        mm2pixels = true;       % flag for conversion of input crop region to pixels

        if ~isa(this,'sctdata')  error('Unknown data type');             end
        if nargin < 2
            [Icrop, origin, cropsize] = crop(this.data);
            x0 = origin(2);     xd = cropsize(2);
            y0 = origin(1);     yd = cropsize(1);
            z0 = origin(3);     zd = cropsize(3);
            mm2pixels = false;
            
        elseif nargin < 7
            error('Missing input arguments. Example call: cropobj = crop(this, x0, y0, z0, xd, yd, zd)');
        end

        % get z crop args
        if isempty(z0),     z0 = 1;                     end
        if isempty(zd),     zd = size(this.data,3);      end
        if zd < 1           zd = 1;                     end

        % get crop args in pixels
        if mm2pixels
            x0 = round(x0/this.voxelsize(1));    xd = round(xd/this.voxelsize(1));
            y0 = round(y0/this.voxelsize(2));    yd = round(yd/this.voxelsize(2));
            z0 = round(z0/this.voxelsize(3));    zd = round(zd/this.voxelsize(3));
        end

        % crop input data
        if size(this.data,3) == 1
            disp('Input data is 2D image; z0 and zd ignored!');
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, 1);
            if isprop(this, 'mask')
                if ~isempty(this.mask)  cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, 1);
                else                 	cropmask = [];
                end
            end
        else
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
            if isprop(this, 'mask')
                if ~isempty(this.mask) 	cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
                else                    cropmask = [];
                end
            end
        end
        
        if nargout > 0
            % create deep copy of input sctdata object (handle is given)
            cropobj = this.copy;
            cropobj.setData(cropdata);
            if ~isempty(cropmask)   cropobj.setMask(cropmask);     end
            cropobj.transformmatrix = [];
            cropobj.centerofrotation = [];            
            if ~isempty(cropobj.offset)
                cropobj.setOffset(this.offset + [x0-1 y0-1 z0-1].*cropobj.voxelsize);
            else
                cropobj.setOffset([x0-1 y0-1 z0-1].*cropobj.voxelsize);
            end
        else
            % write cropped data to input this and update props
            this.setData(cropdata);
            if ~isempty(cropmask)   this.setMask(cropmask);     end
            if ~isempty(this.offset)
                this.setOffset(this.offset + [x0-1 y0-1 z0-1].*this.voxelsize);
            else
                this.setOffset([x0-1 y0-1 z0-1].*this.voxelsize);
            end
        end
        
    end
    
    function [cropobj] = crop_pixel(this, x0, y0, z0, xd, yd, zd)
    %cropobj = crop_pixel(this, x0, y0, z0, xd, yd, zd)
    %   sct crop_pixel fun
    %
    %   WARNING: offset not implemented yet.. origin information will be lost
    %
    %   inputs are the first crop point (x0, y0, z0) and size of the crop region (xd, yd, zd), IN PIXELS
    %
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   pippo.crop_pixel(150, 200, 250, 200, 200, 1).showslice
    %                   pippo.crop_pixel('bbox').showslice
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   29/10/2016
    %   Last update:  09/03/2019
    %   ______________________________________________________

        if ~isa(this,'sctdata')	error('Unknown data type');             end
        if nargin<7 && nargin>2 error('Wrong call.. Example calls:\n\tsctdata.crop_pixel(x0, y0, z0, xd, yd, zd)\n\tsctdata.crop_pixel(''auto'')'); end
        switch nargin
            case 7
                if isempty(x0)  x0 = 1; end
                if isempty(y0)  y0 = 1; end
                if isempty(z0)  z0 = 1; end
                if isempty(xd)  xd = this.size(1)-x0+1;  end
                if isempty(yd)  yd = this.size(2)-y0+1;  end
                if isempty(zd)  zd = this.size(3)-z0+1;  end
                
            case 2
                if isempty(x0)
                    x0 = 'bbox';
                    warning('crop_pixel mode:  bbox');
                end
                validateattributes(x0, {'char'}, {'nonempty'}, mfilename, 'mode', 1);
                x0 = validatestring(x0,{'bbox'});
                if ~isempty(this.mask)
                    [y0, yd, x0, xd, z0, zd] = bbox(this.mask, 1);
                else
                    [y0, yd, x0, xd, z0, zd] = bbox(threshbone(this.data), 1);
                end
        end
        
        %% crop input data
        if size(this.data,3) == 1
            disp('Input data is 2D image; z0 and zd ignored!');
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, 1);
            if isprop(this, 'mask')
                if ~isempty(this.mask)  cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, 1);
                else                 	cropmask = [];
                end
            end
        else
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
            if isprop(this, 'mask')
                if ~isempty(this.mask) 	cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
                else                    cropmask = [];
                end
            end
        end

        if nargout > 0
            %% create deep copy of input sctdata object (handle is given)
            cropobj = this.copy;
            cropobj.setData(cropdata);
            if ~isempty(cropmask)	cropobj.setMask(cropmask);  end
            if ~isempty(cropobj.offset)
                cropobj.setOffset(this.offset + [x0-1 y0-1 z0-1].*cropobj.voxelsize);
            else
                cropobj.setOffset([x0-1 y0-1 z0-1].*cropobj.voxelsize);
            end
            
        else
            %% write cropped data to input this and update props
            this.setData(cropdata);
            if ~isempty(cropmask)	this.setMask(cropmask);     end
            if ~isempty(this.offset)
                this.setOffset(this.offset + [x0-1 y0-1 z0-1].*this.voxelsize);
            else
                this.setOffset([x0-1 y0-1 z0-1].*this.voxelsize);
            end
            
        end

    end
    
    function [cropobj] = crop_pixel_fromcenter(this, xd, yd, zd)
    %cropobj = crop_pixel_fromcenter(this, xd, yd, zd)
    %   crop sctdata from center of image with given crop size
    %   inputs are the size of the crop region (xd, yd, zd), IN PIXELS
    %
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   pippo.crop_pixel(200, 200, 1).showslice
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   27/12/2016
    %   Last update:  27/12/2016
    %   ______________________________________________________
    
        warning('sctdata.crop: offset not implemented yet.. origin information will be lost!');
    
        if ~isa(this,'sctdata')	error('Unknown data type');             end
        if nargin < 4               error('Missing input arguments. Example call: cropobj = crop_pixel_fromcenter(this, xd, yd, zd)');             end
        
        %% get x0, y0, z0
        tmp = round(size(this.data)/2);
        tmp = [tmp(2) tmp(1) tmp(3)]-ceil([xd yd zd]/2)+1;
        x0 = tmp(1);    y0 = tmp(2);    z0 = tmp(3);
        
        %% check bounds
        if z0 <= 0
            fprintf('zd > z_size (%i > %i); z0 = 1\tzd = %i\n', zd, size(this.data,3), size(this.data,3));
            z0 = 1;     zd = size(this.data,3);
        end
        
        if y0 <= 0
            fprintf('yd > y_size (%i > %i); y0 = 1\tyd = %i\n', yd, size(this.data,1), size(this.data,1));
            y0 = 1;     yd = size(this.data,1);
        end
        
        if x0 <= 0
            fprintf('xd > x_size (%i > %i); x0 = 1\txd = %i\n', xd, size(this.data,2), size(this.data,2));
            x0 = 1;     xd = size(this.data,2);
        end

        %% crop input data
        if size(this.data,3) == 1
            disp('Input data is 2D image; zd ignored!');
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, 1);
            if isprop(this, 'mask')
                if ~isempty(this.mask)  cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, 1);
                else                 	cropmask = [];
                end
            end
        else
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
            if isprop(this, 'mask')
                if ~isempty(this.mask) 	cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
                else                    cropmask = [];
                end
            end
        end

        if nargout > 0
            %% create deep copy of input sctdata object (handle is given)
            cropobj = this.copy;
            cropobj.data = cropdata;
            cropobj.mask = cropmask;
            if ~isempty(cropobj.offset)
                cropobj.setOffset(this.offset + [x0-1 y0-1 z0-1].*cropobj.voxelsize);
            else
                cropobj.setOffset([x0-1 y0-1 z0-1].*cropobj.voxelsize);
            end
        else
            %% write cropped data to input this and update props
            this.data = cropdata;
            this.mask = cropmask;
            if ~isempty(this.offset)
                this.setOffset(this.offset + [x0-1 y0-1 z0-1].*this.voxelsize);
            else
                this.setOffset([x0-1 y0-1 z0-1].*this.voxelsize);
            end
        end
    end
    
    function [cropobj] = crop_around_pixel(this, x0, y0, z0, xd, yd, zd)
    %cropobj = crop_around_pixel(this, x0, y0, z0, xd, yd, zd)
    %   crop sctdata around given pixel with given crop size
    %   inputs are the center of the crop region (x0, y0, z0)
    %   and its size (xd, yd, zd), IN PIXELS
    %
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   pippo.crop_pixel(300, 400, 1, 200, 200, 1).showslice
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   27/12/2016
    %   Last update:  27/12/2016
    %   ______________________________________________________
    
        warning('sctdata.crop: offset not implemented yet.. origin information will be lost!');
    
        if ~isa(this,'sctdata')	error('Unknown data type');             end
        if nargin < 7               error('Missing input arguments. Example call: cropobj = crop_around_pixel(this, x0, y0, z0, xd, yd, zd)');             end
        
        if isempty(z0),     z0 = 1;             end
        if isempty(zd),     zd = this.size(3);  end
        
        %% get x0, y0, z0
        tmp = [x0 y0 z0] - ceil([xd yd zd]/2)+1;
        x0_ = tmp(1);    y0_ = tmp(2);    z0_ = tmp(3);
        
        %% check bounds
        if z0_ <= 0
            fprintf('z0-(zd/2) < 0 (%i-(%i/2)=%i); z0 = 1\n', z0, zd, z0_);
            zd = zd + z0_;      z0_ = 1;
        end
        
        if y0_ <= 0
            fprintf('y0-(yd/2) < 0 (%i-(%i/2)=%i); y0 = 1\n', y0, yd, y0_);
            yd = yd + y0_;      y0_ = 1;
        end
        
        if x0_ <= 0
            fprintf('x0-(xd/2) < 0 (%i-(%i/2)=%i); x0 = 1\n', x0, xd, x0_);
            xd = xd + x0_;      x0_ = 1;
        end
        
        if z0_+round(zd/2) > this.size(3)
            fprintf('z0+(zd/2) > this.size(3) (%i>%i); zd = %i\n', z0_+round(zd/2), this.size(3), 2*(this.size(3)-z0_+1));
            zd = 2*(this.size(3)-z0_+1);
        end
        
        if y0_+round(yd/2) > this.size(2)
            fprintf('y0+(yd/2) > this.size(2) (%i>%i); yd = %i\n', y0_+round(yd/2), this.size(2), 2*(this.size(2)-y0_+1));
            yd = 2*(this.size(2)-y0_+1);
        end
        
        if x0_+round(xd/2) > this.size(1)
            fprintf('x0+(xd/2) > this.size(1) (%i>%i); xd = %i\n', x0_+round(xd/2), this.size(1), 2*(this.size(1)-x0_+1));
            xd = 2*(this.size(1)-x0_+1);
        end
        
        x0 = x0_;       y0 = y0_;       z0 = z0_;

        %% crop input data
        if size(this.data,3) == 1
            disp('Input data is 2D image; zd ignored!');
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, 1);
            if isprop(this, 'mask')
                if ~isempty(this.mask)  cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, 1);
                else                 	cropmask = [];
                end
            end
        else
            cropdata = this.data(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
            if isprop(this, 'mask')
                if ~isempty(this.mask) 	cropmask = this.mask(y0:y0+yd-1, x0:x0+xd-1, z0:z0+zd-1);
                else                    cropmask = [];
                end
            end
        end

        if nargout > 0
            %% create deep copy of input sctdata object (handle is given)
            cropobj = this.copy;
            cropobj.data = cropdata;
            cropobj.mask = cropmask;
            if ~isempty(cropobj.offset)
                cropobj.setOffset(this.offset + [x0-1 y0-1 z0-1].*cropobj.voxelsize);
            else
                cropobj.setOffset([x0-1 y0-1 z0-1].*cropobj.voxelsize);
            end
        else
            %% write cropped data to input this and update props
            this.data = cropdata;
            this.mask = cropmask;
            if ~isempty(this.offset)
                this.setOffset(this.offset + [x0-1 y0-1 z0-1].*this.voxelsize);
            else
                this.setOffset([x0-1 y0-1 z0-1].*this.voxelsize);
            end
        end
    end
    
    function [resobj] = resample(this, resamplefactor)
    %resobj = resample(this, resamplefactor)
    %   sct resample method
    %
    %   the resample method doesn't implement interpolation for 3D data
    %   if you are running a MATLAB version >= 2017a you should use imresize3
    %
    %   resamples sctdata.data with given scaling resamplefactor
    %
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   pippo_res = pippo.resample(10)
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   29/10/2016
    %   Last update:  14/11/2016
    %   ______________________________________________________
    
        if ~isa(this,'sctdata')     error('Unknown data type');             end
        if nargin < 2               error('resamplefactor not given. Example call: resobj = resample(resamplefactor)');             end
        if resamplefactor < 1       error('resamplefactor must be greater than 1! RESAMPLE only allows image downsample. for image resize use INTERPOLATE');     end

        %% round resample resamplefactor
        resamplefactor = round(resamplefactor);

        %% initialize output resampled object
        resobj = sctdata;
        resobj.voxelsize                = this.voxelsize*resamplefactor;
        resobj.elementtype              = this.elementtype;
        resobj.rawfile                  = this.rawfile;
        resobj.headerfile               = this.headerfile;
        resobj.transformmatrix          = this.transformmatrix;
        resobj.centerofrotation         = this.centerofrotation;
        resobj.offset                   = this.offset;
        resobj.dataoffset               = this.dataoffset;
        resobj.intensity                = this.intensity;
        resobj.energy                   = this.energy;
        resobj.mu_scaling               = this.mu_scaling;
        resobj.slope_dens               = this.slope_dens;
        resobj.offset_dens              = this.offset_dens;
        resobj.nr_of_projections        = this.nr_of_projections;
        resobj.binarydata               = this.binarydata;
        resobj.thresholdvalue           = this.thresholdvalue;
        resobj.BMDdata                  = this.BMDdata;

        %% resample input data
        if this.size(3) == 1                 % 2D DATA: USE MATLAB IMRESIZE
            resobj.data = imresize(this.data, 1/resamplefactor, 'box');
            resobj.dims = 2;
        else                                % 3D DATA
            
            warning('sctdata.resample: does not implement interpolation for 3D data!');
            
            x_last = floor(size(this.data,2)/resamplefactor);
            y_last = floor(size(this.data,1)/resamplefactor);
            z_last = floor(size(this.data,3)/resamplefactor);

            resobj.data = zeros(y_last, x_last, z_last, class(this.data));
            fprintf('Start resample with factor: %i\n',resamplefactor);
            fprintf('slice:        res');    
            tic
            for kk=1:z_last
                z_start = (kk-1)*resamplefactor+1;
                z_end = kk*resamplefactor;
                data_to_resample = this.data(:, :, z_start:z_end);
                fprintf(1,[repmat('\b',1,numel(num2str(kk))+numel(num2str(z_last))) '\b%i/%i'], kk, z_last);
                for ii=1:x_last
                    x_start = (ii-1)*resamplefactor+1;
                    x_end = ii*resamplefactor;
                    for jj=1:y_last
                        y_start = (jj-1)*resamplefactor+1;
                        y_end = jj*resamplefactor;
                        resobj.data(jj, ii, kk) = mean(mean(mean(data_to_resample(y_start:y_end, x_start:x_end, :))));
                    end
                end
            end
            fprintf('\n')
            toc

            resobj.dims = 3;
        end

    end
        
    function [calibdata] = BMDcalibrate(this)
    %calibdata = BMDcalibrate(this)
    %   calibrates data to BMD
    %       if output argument is given BMDcalibrate returns new object
    %       with BMD calibrated data. If no output argument is given
    %       BMDcalibrate overwrites data field of input sctdata object.
    %
    %       Private dicom tags from Scanco for BMD calibration
    %   	muScaling = double(dicominfo.Private_0029_1000);
    %       slope_dens = double(dicominfo.Private_0029_1004);
    %       offset_dens = double(dicominfo.Private_0029_1005);
    %
    %       Rule used:  BMDROI = ROI.*(slope_dens/muScaling) + offset_dens;
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   31/10/2016
    %   Last update:  09/03/2019
    %   ______________________________________________________

        if ~isa(this,'sctdata')         error('Unknown data type');         	end

        %% check calibration parameters
        if isempty(this.mu_scaling)     error('mu_scaling field is empty');     end
        if isempty(this.slope_dens)     error('slope_dens field is empty');     end
        if isempty(this.offset_dens)	error('offset_dens field is empty');    end

        %% calibration
        if nargout > 0
            %% create deep copy of input sctdata object (handle is given)
            calibdata = this.copy;
            calibdata.data = int16(this.data).*(this.slope_dens/this.mu_scaling) + this.offset_dens;
            calibdata.min_data_value = min(calibdata.data(:));
            calibdata.max_data_value = max(calibdata.data(:));
            calibdata.BMDdata = true;
            calibdata.setOrthoslices;
        else
            %% overwrite calibrated BMD data to input this
            this.data = int16(this.data).*(this.slope_dens/this.mu_scaling) + this.offset_dens;
            this.min_data_value = min(this.data(:));
            this.max_data_value = max(this.data(:));
            this.BMDdata = true;
            this.setOrthoslices;
            fprintf('Data calibrated to BMD!\n');
        end

    end
    
    function [scaleddata] = scaledata(this, slope, offset)
    %scaleddata = scaledata(this)
    %   linearly scale data with given slope and offset
    %       if output argument is given scaledata returns new object
    %       with scaled data. If no output argument is given
    %       scaledata overwrites data field of input sctdata object.
    %
    %   calibdata = (data * slope) + offset
    %   ______________________________________________________

        if ~isa(this,'sctdata')  error('Unknown data type');         	end
        if nargin < 3,              error('Not enough input arguments!');   end

        %% calibration
        if nargout > 0
            %% create deep copy of input sctdata object (handle is given)
            scaleddata = this.copy;
            scaleddata.data = this.data.*slope + offset;
            scaleddata.min_data_value = min(min(min(scaleddata.data)));
            scaleddata.max_data_value = max(max(max(scaleddata.data)));
        else
            %% overwrite calibrated data of input this
            this.data = this.data.*slope + offset;
            this.min_data_value = min(min(min(this.data)));
            this.max_data_value = max(max(max(this.data)));
            fprintf('Data scaled!\n');
        end

    end
    
    function [polydata] = polydata(this, p)
    %[polydata] = polydata(this, coeff)
    %   scale data with linear model:
    %           dataout = p(1)*datain^2 + p(2)*datain + p(3)
    %   or:
    %           dataout = p(1)*datain + p(2)
    %
    %   according to the number of elements in p.
    %   if output argument is given polydata returns new object
    %   with scaled data. If no output argument is given
    %   polydata overwrites data field of input sctdata object.
    %   ______________________________________________________

        if ~isa(this,'sctdata')  error('Unknown data type');         	end
        if nargin < 3,              error('Not enough input arguments!');   end

        %% calibration
        if nargout > 0
            %% create deep copy of input sctdata object (handle is given)
            scaleddata = this.copy;
            scaleddata.data = this.data.*slope + offset;
            scaleddata.min_data_value = min(scaleddata.data(:));
            scaleddata.max_data_value = max(scaleddata.data(:));
        else
            %% overwrite calibrated data of input this
            this.data = this.data.*slope + offset;
            this.min_data_value = min(this.data(:));
            this.max_data_value = max(this.data(:));
            fprintf('Data scaled!\n');
        end

    end
    
    function [outobj] = setdatarange(this, a, b)
    %outobj = setdatarange(this, a, b)
    %   resets data range to given minimum and maximum
    %
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   pippo.showhist
    %                   pippo.setdatarange(-100, 2000);
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   01/11/2016
    %   Last update:  01/11/2016
    %   ______________________________________________________

        if ~isa(this,'sctdata')  error('Unknown data type');             end
        if nargin < 2               error('Data range not given. Example call: out = setdatarange(-200, 3000)');        end

        %% check input values
        MIN = min([a b]);       MAX = max([a b]);

        %% apply range to data
        data = this.data;
        data(data < MIN) = MIN;
        data(data > MAX) = MAX;

        %% 
        if nargout > 0
            %% create deep copy of input sctdata object (handle is given)
            outobj = this.copy;
            outobj.data = data;         clear data;
            outobj.min_data_value = MIN;
            outobj.max_data_value = MAX;
        else
            %% overwrite input object data field
            this.data = data;            clear data
            this.min_data_value = MIN;
            this.max_data_value = MAX;
        end

    end
    
    function [N, bins] = hist(this, ROI, nbins)
    %[N, bins] = hist(this, ROI, nbins)
    %   sct histogram
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   [N, bins] = pippo.hist              % histogram of whole data
    %                   [N, bins] = pippo.hist(ROI, 256)	% histogram of given ROI
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   09/11/2016
    %   Last update:  09/11/2016
    %   ______________________________________________________

        if nargin < 3,              nbins = 100;                            end
        if ~isa(this,'sctdata')  error('Unknown data type');             end
        if ~isprop(this, 'data')     error('object contains no data');       end

        if nargin<2
            if verLessThan('matlab','8.4')
                % -- Code to run in MATLAB R2014a and earlier here --
                [N, bins] = hist(double(reshape(this.data, [numel(this.data) 1])), nbins);
            else
                % -- Code to run in MATLAB R2014b and later here --
                [N, edges] = histcounts(this.data, nbins);
                bins = edges(1:end-1)+(edges(2)-edges(1));
            end
        else
            if verLessThan('matlab','8.4')
                % -- Code to run in MATLAB R2014a and earlier here --
                [N, bins] = hist(double(reshape(ROI, [numel(ROI) 1])), nbins);
            else
                % -- Code to run in MATLAB R2014b and later here --
                [N, edges] = histcounts(ROI, nbins);
                bins = edges(1:end-1)+(edges(2)-edges(1));
            end
        end

    end
            
    function [new] = copy(this)
    % Make a copy of a handle object.

        % Instantiate new object of the same class.
        new = feval(class(this));

        % Copy all non-hidden properties.
        p = properties(this);
        for i = 1:length(p)
            new.(p{i}) = this.(p{i});
        end
    end
    
    %     function [resobj] = resize(this, sizefactor)
    %     %[resobj] = resize(this, sizefactor)
    %     %   sct resize method
    %     %   resizes data and (if exists) mask of sctdata object with given size factor
    %     %   sizefactor      can be either a scalar or [row col slice] size of
    %     %                   data after resize.
    %     %
    %     %   example call:   pippo = sctdata.load('./pippo.mhd');    % load 2D image
    %     %                   pippo_res = pippo.resize(2)                 % resize with factor 2
    %     %                   pippo_res = pippo.resize([1200 1000])       % resize to new image size
    %    
    %         if ~isa(this,'sctdata') error('sctdata.resize: Unknown data type');	end
    %         if nargin < 2               error('sctdata.resize: sizefactor not given. Example call: resobj = this.resize(sizefactor)');     end
    % %         if sizefactor < 1       error('resamplefactor must be greater than 1! RESAMPLE only allows image downsample. for image resize use INTERPOLATE');     end
    %         
    %         if length(sizefactor) == 2
    %             if this.dims ~= 2       error('sctdata.resize: data must be 2D');	end
    %             
    %             % resize 2D data and mask
    %             size_pre = this.size;
    %             interpdata = imresize(this.data, sizefactor);
    %             if ~isempty(this.mask)
    %                 interpmask = imresize(this.mask, sizefactor);
    %             else
    %                 interpmask = [];
    %             end
    %             voxelsize = this.voxelsize./([sizefactor([2 1])./size_pre(1:2) 1]);
    %         
    %         elseif length(sizefactor) == 3
    %             if this.dims ~= 3       error('sctdata.resize: data must be 3D');	end
    %             if sizefactor(3) == this.size(3)
    %                 % imresize slicewise
    %                 size_pre = this.size;
    %                 
    %                 if ~isempty(this.mask)
    %                     interpdata = sct.imresize(this.data, sizefactor);
    %                     interpmask = sct.imresize(this.mask, sizefactor);
    %                 else
    %                     interpdata = sct.imresize(this.data, sizefactor);
    %                     interpmask = [];
    %                 end
    %                 
    %                 voxelsize = this.voxelsize./(sizefactor([2 1 3])./size_pre);
    %                 
    %             else
    %                 % interp3
    %                 error('sct.sctdata.resize: 3D interpolate not implemented yet');
    %             end
    %         end
    %         
    %         if nargout > 0
    %             %% create deep copy of input sctdata object (handle is given)
    %             resobj = this.copy;
    %             resobj.data = interpdata;
    %             resobj.mask = interpmask;
    %             resobj.voxelsize = voxelsize;
    %             % resobj will lose transformmatrix, centerofrotation and offset of the original data
    %             resobj.transformmatrix = [];
    %             resobj.centerofrotation = [];
    %             resobj.offset = [];
    %         else
    %             %% write resized data to input obj and update props
    %             this.data = interpdata;
    %             this.mask = interpmask;
    %             this.voxelsize = voxelsize;
    %             % this will lose transformmatrix, centerofrotation and offset
    %             this.transformmatrix = [];
    %             this.centerofrotation = [];
    %             this.offset = [];
    %         end
    %     
    %     end

    %     function [bindata] = threshold(this, thresh)
    %     %bindata = threshold(this, thresh)
    %     %   sctdata threshold method
    %     %       if output argument is given threshold returns sctdata
    %     %       object of binary data. If no output argument is given
    %     %       threshold adds binary mask field to input sctdata object.
    %     %
    %     %   the input thresh can be numeric (1 value for single global threshold, array of 2 values for 2-level global threshold)
    %     %   or one of the following string options:
    %     %       'otsu'      for thresholding with the Otsu's method
    %     %       'adaptive'  for thresholding of SAM images with Adaptive threshold Method from Lakshman 2007
    %     %
    %     %   by default (no argument given) threshold uses the Otsu's method
    %     %
    %     %   example call:   pippo = sct.sctdata.load('./pippo.mhd');
    %     %                   pippobin = pippo.threshold                  % threshold with Otsu's method
    %     %                   pippobin = pippo.threshold(128);            % threshold with global single level threshold of 128
    %     %                   pippobin = pippo.threshold([5 12]);         % threshold with global threshold (MIN and MAX)
    %     %                   pippobin = pippo.threshold('adaptive')      % threshold with Adaptive threshold (use for SAM images)
    %     %                   pippo.threshold                             % adds mask obtained with Otsu's method
    %     %                   pippo.threshold.showslice               	% display thresholded mid slice
    %     %   ______________________________________________________
    %     %
    %     %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %     %   BSRT - Charite Berlin
    %     %   Created on:   25/10/2016
    %     %   Last update:  01/02/2017
    %     %   ______________________________________________________
    % 
    %         if ~isa(this,'sct.sctdata')     error('Unknown data type');             end
    %         if nargin < 2                   thresh = [];                            end
    % 
    %         %% initialize binary data matrix
    %         data = false(size(this.data));
    % 
    %         %% apply threshold based on input
    %         if isempty(thresh)
    %             fprintf('No threshold level given.. I will use Otsu method.\n');
    %             thresh = 'otsu';
    %         end
    %         
    %         if isnumeric(thresh)
    %             % global threshold (single or double) given by user
    %             if length(thresh) > 2     error('too many input thresh arguments');
    %             
    %             elseif length(thresh) == 2
    %                 % MIN and MAX global threshold
    %                 fprintf('Thresholding with 2-level global thresh:\n\tmin: %2.2f\n\tmin: %2.2f\n', min(thresh), max(thresh));
    %                 mask1 = this.data > min(thresh);
    %                 mask2 = this.data < max(thresh);
    %                 data(mask1 & mask2) = 1;
    %                 data = ~data;
    %                 
    %             elseif length(thresh) == 1
    %                 % single valued global threshold
    %                 fprintf('Thresholding with single level global thresh:\n\tthresh: %2.2f\n', thresh);
    %                 data(this.data > thresh) = 1; 
    %             end
    %             
    %             threshval = thresh;
    %             
    %         else
    %             switch(thresh)
    %                 case('otsu')
    %                     % Otsu's method
    %                     threshval = multithresh(this.data);                     % obtain threshold level with Otsu's Method
    %                     fprintf('Thresholding with Otsu Method:\n\tthresh: %2.2f\n', threshval);
    %                     data(this.data > threshval) = 1;
    %                     
    %                 case('adaptive')
    %                     threshval = sct.Thresholder3D.adaptivethresh(this.data);
    %                     fprintf('Thresholding with Adaptive threshold Method from Lakshman 2007:\n\tthresh: %2.2f\n', threshval);
    %                     data(this.data > threshval) = 1;
    %                     
    %                 case('iterative')
    %                     threshval = sct.Thresholder3D.isodata(this.data);
    %                     fprintf('Thresholding with iterative isodata method from Ridler et al:\n\tthresh: %2.2f\n', threshval);
    %                     data(this.data > threshval) = 1;
    %                     
    %             end
    %         end
    %         
    %         %% output 
    %         if nargout > 0
    %             %% create deep copy of input sctdata object (handle is given)
    %             bindata = this.copy;
    %             bindata.binarydata = true;
    %             bindata.data = data;
    %             clear data;
    %         else
    %             %% write binary mask to input this
    %             this.mask = data;
    %             this.thresholdvalue = threshval;
    %             clear data
    %         end
    % 
    %     end
    
    %% sctdata visualization methods
    function showslice(this, sln, dir)
    %visualize one slice of volume data
    %   example call:   showslice           % displays Z midslice
    %                   showslice(120)      % displays Z slice 120
    %                   showslice(120,2)    % displays Y slice 120
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   25/10/2016
    %   Last update:  23/01/2018
    %   ______________________________________________________

        if ~isa(this,'sctdata')     error('Unknown data type');             end
        if isempty(this.data)       error('data field is empty');           end

        % if dir input not given displays a Z slice
        if nargin<3                 dir = 3;                                end
        % if sln input not given displays the midslice
        if nargin<2                 sln = round(this.size(dir)/2);          end
        if isempty(sln)             sln = round(this.size(dir)/2);          end

        figure;
        sliceim = this.getSlice(sln, dir);                  % get the slice image
        
        switch dir
            case 1
                if ~isempty(this.voxelsize)
                    % if voxelsize field is present displays image with physical size
                    if length(this.voxelsize) > 1
                        imagesc([1:this.size(2)]*this.voxelsize(2), [1:this.size(3)]*this.voxelsize(3), sliceim);
                    else        % isotropic voxelsize
                        imagesc([1:this.size(2)]*this.voxelsize, [1:this.size(3)]*this.voxelsize, sliceim);
                    end
                    xlabel('y [mm]');
                    ylabel('z [mm]');
                else
                    imagesc(sliceim);
                    xlabel('y [pixels]');
                    ylabel('z [pixels]');
                end
            case 2
                if ~isempty(this.voxelsize)
                    % if voxelsize field is present displays image with physical size
                    if length(this.voxelsize) > 1
                        imagesc([1:this.size(1)]*this.voxelsize(1), [1:this.size(3)]*this.voxelsize(3), sliceim);
                    else        % isotropic voxelsize
                        imagesc([1:this.size(1)]*this.voxelsize, [1:this.size(3)]*this.voxelsize, sliceim);
                    end
                    xlabel('x [mm]');
                    ylabel('z [mm]');
                else
                    imagesc(sliceim);
                    xlabel('x [pixels]');
                    ylabel('z [pixels]');
                end
            case 3
                if ~isempty(this.voxelsize)
                    % if voxelsize field is present displays image with physical size
                    if length(this.voxelsize) > 1
                        imagesc([1:this.size(1)]*this.voxelsize(1), [1:this.size(2)]*this.voxelsize(2), sliceim);
                    else        % isotropic voxelsize
                        imagesc([1:this.size(1)]*this.voxelsize, [1:this.size(2)]*this.voxelsize, sliceim);
                    end
                    xlabel('x [mm]');
                    ylabel('y [mm]');
                else
                    imagesc(sliceim);
                    xlabel('x [pixels]');
                    ylabel('y [pixels]');
                end
        end
        
        colormap gray;
        axis image;
        set(gca,'ydir','normal');

    end
    
    function showslices(this, YZslicen, XZslicen, XYslicen)
        switch nargin
            case 1
                if isempty(this.XYslice)    this.setOrthoslices;    end
                
            case 4
                if isempty(YZslicen), YZslicen = this.YZslicen;     end
                if isempty(XZslicen), XZslicen = this.XZslicen;     end
                if isempty(XYslicen), XYslicen = this.XYslicen;     end
                
                this.setXYslicen(XYslicen);
                this.setXZslicen(XZslicen);
                this.setYZslicen(YZslicen);
                this.setOrthoslices;
        end
        
        figure('units','normalized','outerposition',[0 0 1 1]);
        green = cat(3, zeros(size(this.XYslicemask)), ones(size(this.XYslicemask)), zeros(size(this.XYslicemask)));
        subplot 131;    h=imagesc(this.XYslice);  colormap gray;  axis image;   hold on;    h = imagesc(green);     hold off;   set(h, 'AlphaData', 0.3*double(this.XYslicemask));
        green = cat(3, zeros(size(this.XZslicemask)), ones(size(this.XZslicemask)), zeros(size(this.XZslicemask)));
        subplot 132;    h=imagesc(this.XZslice);  colormap gray;  axis image;   hold on;    h = imagesc(green);     hold off;   set(h, 'AlphaData', 0.3*double(this.XZslicemask));
        green = cat(3, zeros(size(this.YZslicemask)), ones(size(this.YZslicemask)), zeros(size(this.YZslicemask)));
        subplot 133;    h=imagesc(this.YZslice);  colormap gray;  axis image;   hold on;    h = imagesc(green);     hold off;   set(h, 'AlphaData', 0.3*double(this.YZslicemask));

    end
    
    function showmaskslice(this, slicenum)
    %showmaskslice(this, slicenum)
    %   sct binary mask slice visualizator
    %   example call:   pippo = sctdata.load('./pippo.mhd');
    %                   pippo.threshold;
    %                   pippo.showmaskslice         % displays midslice
    %                   pippo.showmaskslice(120)	% displays slice 120
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   02/11/2016
    %   Last update:  02/11/2016
    %   ______________________________________________________

        if ~isa(this,'sctdata')     error('Unknown data type');             end
        if isempty(this.mask)       error('binary mask is empty');          end

        if nargin<2                 slicenum = round(size(this.data,3)/2);   end

        figure;
        if ~isempty(this.voxelsize)
            if length(this.voxelsize) > 1
                imagesc([1:size(this.data,2)]*this.voxelsize(1), [1:size(this.data,1)]*this.voxelsize(2), this.mask(:,:,slicenum));
            else
                imagesc([1:size(this.data,2)]*this.voxelsize(1), [1:size(this.data,1)]*this.voxelsize(2), this.mask(:,:,slicenum));
            end
            xlabel('x [mm]');
            ylabel('y [mm]');
        else
            imagesc(this.mask(:,:,slicenum));
            xlabel('x [pixels]');
            ylabel('y [pixels]');
        end
        colormap gray;
        axis image;
        set(gca,'ydir','normal');

    end
    
    function showhist(this, ROI, nbins)
        %showhist(this, ROI, nbins)
        %   sct histogram visualizator
        %   example call:   pippo = sctdata.load('./pippo.mhd');
        %                   pippo.showhist              % histogram of whole data
        %                   pippo.showslice(ROI, 30)	% histogram of given ROI
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianluca.iori@charite.de)
        %   BSRT - Charite Berlin
        %   Created on:   01/11/2016
        %   Last update:  09/11/2016
        %   ______________________________________________________

        if nargin < 3,              nbins = 256;                            end
        if ~isa(this,'sctdata')  error('Unknown data type');             end
        if ~isprop(this, 'data')     error('object contains no data');       end

        if nargin>=2 && ~isempty(ROI)
            if verLessThan('matlab','8.4')
                % -- Code to run in MATLAB R2014a and earlier here --
                figure; hist(single(reshape(ROI, [numel(ROI) 1])), nbins);
            else
                % -- Code to run in MATLAB R2014b and later here --
                figure; histogram(ROI, nbins);
            end
            xlabel('a.u.');
        else
            if verLessThan('matlab','8.4')
                % -- Code to run in MATLAB R2014a and earlier here --
                figure; hist(single(reshape(this.data, [numel(this.data) 1])), nbins);
            else
                % -- Code to run in MATLAB R2014b and later here --
                figure; histogram(this.data, nbins);
            end
            if this.BMDdata
                xlabel('BMD [mg_{HA}/cc]');
            else
                xlabel('a.u.');
            end
        end
        ylabel('counts');

    end
    
    %% input/output
    function load(this, filename)
        % function load(this, filename)
        % load sctdata
        %
        % supported file formats:    - DICOM data           (*.dcm,*.DCM)
        %                            - MetaImage data       (*.mhd,*.MHD)
        %                                                   (*.mha,*.MHA)
        %                            - (Scanco) ISQ data    (*.isq,*.ISQ)
        %                            - (Scanco) AIM data    (*.aim,*.AIM)
        %                            - VolFloat32 data      (*.vol)
        %                            - MATLAB data          (*.mat,*.MAT)
        %                            - image stacks         (*.jpg;*.gif;*.png;*.bmp;*.tif)
        %
        % example call:     pippo = sct.sctdata
        %                   pippo.load

       % get input file if no argument is given
       if nargin < 2

           [filename, pathname] = uigetfile({  '*.dcm;*.DCM;*.mhd;*.MHD;*.isq;*.ISQ;*.aim;*.AIM;*.mat;*.MAT;*.jpg;*.gif;*.png;*.bmp;*.tif','Supported Files (*.dcm,*.DCM,*.mhd,*.MHD,*.ISQ,*.mat,*.MAT,*.jpg,*.gif,*.png,*.bmp,*.tif)';...
                                               '*.dcm;*.DCM','DICOM Files (*.dcm,*.DCM)';...
                                               '*.mhd;*.MHD','MetaImage data (*.mhd,*.MHD)';...
                                               '*.isq;*.ISQ','Scanco ISQ data (*.isq,*.ISQ)';...
                                               '*.aim;*.AIM','Scanco AIM data (*.aim,*.AIM)';...
                                               '*.vol;*.VOL','Synchrotron CT Data Files (*.vol,*.VOL)';...
                                               '*.mat;*.MAT','MATLAB data (*.mat,*.MAT)';...
                                               '*.jpg;*.gif;*.png;*.bmp;*.tif','Image Files';...
                                               '*.*', 'All Files (*.*)'},'Load sctdata image','MultiSelect','on');
           if isequal(filename,0) || isequal(pathname,0)
                % user pressed cancel
                return;
           end

           if isempty(filename)
               if length(filename) > 1
                   if ~ischar(filename{1})
                       this = [];
                       fprintf('opening abborted');
                       msg  = 'opening abborted';
                       return;
                   end
               elseif ~ischar(filename)
                   this = [];
                   fprintf('opening abborted');
                   msg  = 'opening abborted';
                   return;
               end
           end 
           % multiple file select
           if iscell(filename)
               for i=1:length(filename)
                   filenamecell{i} = [pathname filename{i}];
               end
           else
               filenamecell{1} = [pathname filename];
           end
       else
           % multiple file select
           if iscell(filename)
               for i=1:length(filename)
                    filenamecell{i} = filename{i};
               end
           else
               filenamecell{1} = filename;
           end
       end

        % get input fileparts
        [pathstr,name,ext] = fileparts(filenamecell{1});

        % launch specific load method based on file extension
        switch ext
            case {'.mhd','.MHD'}
                if iscell(filename) && length(filename) > 1
                    this.loadMHD(filenamecell);
                else
                    this.loadMHD(filenamecell{1});
                end
            case {'.isq','.ISQ'}
                this.loadISQ(filenamecell{1});
            case {'.aim','.AIM'}
                this.loadAIM(filenamecell{1});
            case {'.dcm','.DCM'}
                this.loadDICOM(filenamecell);
            case {'.vol','.VOL'}
                this.loadVolFloat32(filenamecell);
            case {'.mat','.MAT'}
                this.loadMAT(filenamecell{1});
            case {'.jpg','.gif','.png','.bmp','.tif'}
                this.loadstack(filenamecell);
        end
    end
    
    function loadMHD(this, filename, slicerange)
        %sct loadMHD fun
        %   slicerange allows to specify the slices to be read
        %
        %   example call:   pippo =  sctdata
        %                   pippo.loadMHD('./pippo.mhd');
        %                   pippo.loadMHD('./pippo.mhd',[20:30]);
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianluca.iori@charite.de)
        %   BSRT - Charite Berlin
        %   Created on:   13/10/2016
        %   Last update:  20/07/2017
        %   ______________________________________________________

        % load MHD and set data
        if nargin < 3,
            slicerange = [];
            if nargin < 2
                [V, header, filename] = MHDdata.load;
            else
                [V, header] = MHDdata.load(filename(1:end-4));
            end
        else
            [V, header] = MHDdata.load(filename, slicerange);
        end
        
        % set data
        if ~isempty(V)      this.setData(V);        end

        %% load mhd header and set header info
        meta_var                = textscan(header, '%s%s', 'delimiter', '=');                               % read header file
        this.setDims(str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'NDims'))))));
        this.setSize(str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'DimSize'))))));
        this.voxelsize          = str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'ElementSpacing')))));
        this.elementtype        = char(meta_var{2}(find(strcmp(deblank(meta_var{1}),'ElementType'))));
        this.rawfile            = char(meta_var{2}(find(strcmp(deblank(meta_var{1}),'ElementDataFile'))));
        this.headerfile         = filename;
        if ~isempty(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'TransformMatrix')))))
            this.setTransformmatrix(str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'TransformMatrix'))))));
        end
        if ~isempty(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'CenterOfRotation')))))
            this.centerofrotation   = str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'CenterOfRotation')))));
        end
        this.offset             = str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'Offset')))));

    end

    function loadDICOM(this, filename)
    %sctdata loadDICOM fun
    %   example call:   pippo =  sctdata
    %                   pippo.loadDICOM;
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   31/10/2016
    %   Last update:  31/01/2018
    %   ______________________________________________________

        if nargin > 1                                                      % input files given
            if isstr(filename)      filename = cellstr(filename);	end    % create cell array
            [data, info] = DICOMdata.load(filename);                       % load given input files
        else                                                               % no input file
            [data, info] = DICOMdata.load;                                 % load DICOMs
        end

        this.setData(data);             % set data
        
        % set header info
        if info.Slices == 1     this.dims = 2;
        elseif info.Slices > 1  this.dims = 3;
        end
        this.size                = [info.Columns info.Rows info.Slices];
        this.voxelsize           = [info.PixelSpacing' info.SliceThickness];
        this.elementtype         = class(this.data);
        this.rawfile             = info.Filename;
        this.headerfile          = info.Filename;
        if isfield(info.PatientName, 'FamilyName')  this.samplename = info.PatientName.FamilyName;  end
        % this.intensity           = headerinfo.intensity;
        % this.energy              = headerinfo.energy;
        if isfield(info,'Private_0029_1000')        this.mu_scaling = double(info.Private_0029_1000);   end
        if isfield(info,'Private_0029_1004')        this.slope_dens = double(info.Private_0029_1004);   end
        if isfield(info,'Private_0029_1005')        this.offset_dens = double(info.Private_0029_1005);  end
        
        % this.nr_of_projections	= headerinfo.nr_of_projections;
        this.min_data_value      = min(min(min(this.data)));
        this.max_data_value      = max(max(max(this.data)));
        this.nr_of_bytes         = info.FileSize;

    end
    
    function loadISQ(this, filename, varargin)
        %function loadISQ(this, filename, x_min, y_min, z_min, size_x, size_y, size_z)
        %subvol = readISQ(this, filename, x_min, y_min, z_min, size_x, size_y, size_z)
        %   sctdata loadISQ method
        %   the method is intended for reading sub portion of large ISQ volume files
        %   the subvolume to be read must be defined through the input parameters:
        %   x_min, y_min, z_min     coordinates of first pixel
        %   size_x, size_y, size_z  size (in pixels) to be read
        %
        %   example call:   pippo =  sctdata
        %                   pippo.readISQ('\\charite.de\centren\#Charite-Central\BCRT\AG-raum-qbam-archiv-read\2015.003.qbam.TaCoSound\data\XtremeCT-II\femur\1955_L\C0001577.ISQ');
        %                   
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianthk.iori@gmail.com)
        %   BSRT - Charite Berlin
        %   Created on:   22/10/2016
        %   Last update:  10/01/2018
        %   ______________________________________________________
    
        % Check the number of input arguments.
        narginchk(1, 9);

        switch nargin
            case 1          % no input
                [data, header, filename] = ISQdata.load;
                this.setData(data);
                
            case 2          % filename
                [data, header, filename] = ISQdata.load(filename);
                this.setData(data);

            case 3          % filename and read mode
                readmode = varargin{1};
                validateattributes(readmode,   {'logical','numeric'}, {'binary'}, mfilename,'readmode', 2);
                if readmode
                    [data, header, filename] = ISQdata.load(filename);
                    this.setData(data);

                else
                    % fetch ISQ header info
                    header = ISQdata.readheader(filename);
                end

            case 8          % filename and x,y,z read bounds
                x_min = varargin{1};
                y_min = varargin{2};
                z_min = varargin{3};
                size_x = varargin{4};
                size_y = varargin{5};
                size_z = varargin{6};
                
                [data, header, filename] = ISQdata.load(filename, x_min, y_min, z_min, size_x, size_y, size_z);
                this.setData(data);

            case 9          % filename, read mode and x,y,z read bounds
                readmode = varargin{1};
                validateattributes(readmode,   {'logical','numeric'}, {'binary'}, mfilename,'readmode', 2);
                if readmode
                    x_min = varargin{2};
                    y_min = varargin{3};
                    z_min = varargin{4};
                    size_x = varargin{5};
                    size_y = varargin{6};
                    size_z = varargin{7};
                    [data, header, filename] = ISQdata.load(filename, x_min, y_min, z_min, size_x, size_y, size_z);
                    this.setData(data);

                else
                    warning('loadISQ: readmode was set to ''FALSE'', read bounds will be ignored!');
                    % fetch ISQ header info
                    header = ISQdata.readheader(filename);
                end

            otherwise
                warning('loadISQ: Input not supported.. Example call: loadISQ(filename, mode, x_min, y_min, z_min, size_x, size_y, size_z)');
                
        end
        
        % set ISQ header info
%         this.setSize([header.x_dim header.y_dim header.z_dim]);
        this.setVoxelsize(1e-3*[header.x_dim_um/header.x_dim header.y_dim_um/header.y_dim header.z_dim_um/header.z_dim]);
        this.elementtype        = header.type;
        this.rawfile            = filename;
        this.headerfile         = filename;
        this.samplename         = header.samplename;
        this.dataoffset         = header.offset;
        this.intensity          = header.intensity;
        this.energy             = header.energy;
        this.mu_scaling         = header.mu_scaling;
        this.nr_of_projections  = header.nr_of_projections;
        this.min_data_value     = header.min_data_value;
        this.max_data_value     = header.max_data_value;
        this.nr_of_bytes        = header.nr_of_bytes;
       
    end
    
    function loadVolFloat32(this, filename)
    %sctdata loadVolFloat32 fun
    %   example call:   pippo =  sctdata
    %                   pippo.loadVolFloat32;
    %   ______________________________________________________
    %
    %   Author: Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   19/01/2018
    %   Last update:  19/01/2018
    %   ______________________________________________________

        % load and set data
        if nargin < 2
            [V, header, filename] = VOLdata.load;
        else
            [V, header] = VOLdata.load(filename);
        end
        
        this.setData(V);        % set data

        %% set header info
        meta_var                = textscan(header, '%s%s', 'delimiter', '=');                               % read header file
        this.setDims(str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'NDims'))))));
        this.setSize(str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'DimSize'))))));
        this.voxelsize          = str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'ElementSpacing')))));
        this.elementtype        = char(meta_var{2}(find(strcmp(deblank(meta_var{1}),'ElementType'))));
        this.rawfile            = char(meta_var{2}(find(strcmp(deblank(meta_var{1}),'ElementDataFile'))));
        this.headerfile         = filename;
        this.transformmatrix    = str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'TransformMatrix')))));
        this.centerofrotation   = str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'CenterOfRotation')))));
        this.offset             = str2num(cell2mat(meta_var{2}(find(strcmp(deblank(meta_var{1}),'Offset')))));

    end
    
    function resampleISQ(this, filename, resfactor)
        %RESAMPLEISQ resample and load ISQ data
        %   the method is intended for reading and downsampling of large ISQ volume files
        %   
        %   example call:   pippo =  sctdata
        %                   pippo.resampleISQ('\\charite.de\centren\#Charite-Central\BCRT\AG-raum-qbam-archiv-read\2015.003.qbam.TaCoSound\data\XtremeCT-II\femur\1955_L\C0001577.ISQ', 10);
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianluca.iori@charite.de)
        %   BSRT - Charite Berlin
        %   Created on:   21/11/2016
        %   Last update:  10/01/2018
        %   ______________________________________________________

        if nargin < 2,   disp('resampleISQ: No input file..');    return;   end

        % fetch ISQ header info
        header = ISQdata.readheader(filename);
        
        % set header info
        % this.size                = [header.x_dim header.y_dim header.z_dim];
        this.setVoxelsize(1e-3*[header.x_dim_um/header.x_dim header.y_dim_um/header.y_dim header.z_dim_um/header.z_dim]);
        this.elementtype        = header.type;
        this.rawfile            = filename;
        this.headerfile         = filename;
        this.samplename         = header.samplename;
        this.dataoffset         = header.offset;
        this.intensity          = header.intensity;
        this.energy             = header.energy;
        this.mu_scaling         = header.mu_scaling;
        this.nr_of_projections  = header.nr_of_projections;
        this.min_data_value     = header.min_data_value;
        this.max_data_value     = header.max_data_value;
        this.nr_of_bytes        = header.nr_of_bytes;

        % read and downsample ISQ data
        this.setData(ISQdata.resample(filename, resfactor));
        
    end
    
    function loadAIM(this, filename)
    
        if nargin < 2,   disp('loadAIM: No input file..');    return;   end

        [header, data] = AIMdata.load(filename);

        this.setData(data);
%         this.voxelsize           = 1e-3*[headerinfo.x_dim_um/headerinfo.x_dim headerinfo.y_dim_um/headerinfo.y_dim headerinfo.z_dim_um/headerinfo.z_dim];
        this.elementtype         = header.type;
        this.rawfile             = filename;
        this.headerfile          = filename;
%         this.samplename          = header.samplename;
        this.dataoffset          = header.offset;
%         this.intensity           = header.intensity;
%         this.energy              = header.energy;
%         this.mu_scaling          = header.mu_scaling;
%         this.nr_of_projections	= header.nr_of_projections;
%         this.min_data_value      = header.min_data_value;
%         this.max_data_value      = header.max_data_value;
%         this.nr_of_bytes         = header.nr_of_bytes;

    end
      
    function loadstack(this, filename)
    %loadstack(this, filename)
    %   sct read image sequence method
    %   example call:   pippo =  sctdata
    %                   pippo.loadstack;
    %   ______________________________________________________
    %
    %   Created on:   04/11/2016
    %   Last update:  10/01/2018
    %   ______________________________________________________

        if nargin > 1                                                           % input files given
            if isstr(filename)      filename{1} = filename;      end            % create cell array
            [data, info, filename] = imagestackdata.load(filename);             % load given input files
        else                                                                    % no input file
            [data, info, filename] = imagestackdata.load;                       % load image sequence
        end

        % set data
        this.setData(data);
        
        % set header info
        if info.Slices == 1     this.dims = 2;
        elseif info.Slices > 1  this.dims = 3;
        end
        %         this.size                = [info.Columns info.Rows info.Slices];
        if isfield(info,'ImagePixelSizeum')
            this.voxelsize           = [info.ImagePixelSizeum info.ImagePixelSizeum info.ImagePixelSizeum]*1e-3;
        else
            this.voxelsize           = [];
        end
        this.elementtype         = info.elementtype;
        this.rawfile             = filename{1};
        this.headerfile          = filename{1};
        if isfield(info,'ImagePixelSizeum')
            this.intensity           = info.SourceCurrentuA;
        else
            this.intensity           = [];
        end
        if isfield(info,'ImagePixelSizeum')
            this.energy              = info.SourceVoltagekV;
        else
            this.energy              = [];
        end
        this.min_data_value      = min(min(min(this.data)));
        this.max_data_value      = max(max(max(this.data)));
        % if isfield(info.PatientName, 'FamilyName')  this.samplename = info.PatientName.FamilyName;  end
        % this.mu_scaling          = double(info.Private_0029_1000);
        % this.slope_dens          = double(info.Private_0029_1004);
        % this.offset_dens         = double(info.Private_0029_1005);
        % this.nr_of_projections	= headerinfo.nr_of_projections;
        % this.nr_of_bytes         = info.FileSize;

    end
        
    function readWAVE(this, filename)
    %function readWAVE(this, filename)
    %sct read WAVE matlab file
    %   example call:   pippo =  sctdata
    %                   pippo.readWAVE('./pippo.mat');
    %   ______________________________________________________
    %
    %   Author:     Gianluca Iori (gianluca.iori@charite.de)
    %   BSRT - Charite Berlin
    %   Created on:   02/11/2016
    %   Last update:  19/12/2016
    %   ______________________________________________________

        %% get filename if no file is given as input
        if nargin < 2
            [filename, pathname] = uigetfile({'*.mat;*.MAT','MATLAB data (*.mat,*.MAT)';'*.*', 'All Files (*.*)'},'Load MATLAB WAVE data');
            filename = [pathname filename];
        end

        %% load data
        onda = load(filename);
        
        %% set header info and data to object
        if isfield(onda,'Result')
            % Bonhoff data structure format
            this.data = flipud(onda.Result.Z);
            this.dims = 2;
            this.size = [size(onda.Result.Z,2) size(onda.Result.Z,1) 1];
            this.voxelsize = [onda.Result.X(2)-onda.Result.X(1) onda.Result.Y(2)-onda.Result.Y(1) 1];
            this.elementtype = char(class(onda.Result.Z));
            this.rawfile = filename;
        elseif isfield(onda,'WAVE')
            % Typical NewSAM WAVE data structure format
            this.rawfile = filename;
            this.dims = 2;
            this.voxelsize = [onda.WAVE.X(2)-onda.WAVE.X(1) onda.WAVE.Y(2)-onda.WAVE.Y(1) 1];
            if ~isempty(onda.WAVE.IMP)
                this.data = flipud(onda.WAVE.IMP);
                this.size = [size(onda.WAVE.IMP,2) size(onda.WAVE.IMP,1) 1];
                this.elementtype = char(class(onda.WAVE.IMP));
            else
                warning('IMP field is empty!');
            end
        else
            error('unknown MATLAB data structure');
        end
     
    end
    
    function save(this, filename)
    % function save(this, filename)
    % save sctdata
    %
    % supported file formats:    - MetaImage data    (*.mhd,*.MHD)

       % get output filename if no argument is given
       if nargin < 2
            [filename, pathname] = uiputfile({ '*.mhd;*.MHD;*.mat;*.MAT;*.png;*.PNG','Supported Files (*.mhd,*.MHD,*.mat,*.MAT,*.png,*.PNG)';...
                                               '*.mhd;*.MHD','MetaImage data (*.mhd,*.MHD)';...
                                               '*.mat;*.MAT','MATLAB data (*.mat,*.MAT)';...
                                               '*.png;*.PNG','.png MID-slices';...
                                               '*.*', 'All Files (*.*)'},'Save sctdata image');
            filename = [pathname filename];
       end
       
       % get input fileparts
       [pathstr,name,ext] = fileparts(filename);
       
       % launch spcific save method
        switch ext
            case {'.mhd','.MHD'}
                this.writeMHD(filename);
            case {'.png','.PNG'}
                this.writeMidplanes(filename);
            case {'.mat','.MAT'}
                this.writeMAT(filename);
        end
       
    end
    
    function writeMHD(this, filename)
        %function writeMHD(this, filename)
        %   sct write MHD fun
        %   example call:   pippo.writeMHD('./pippo.mhd');
        %   ______________________________________________________
        %
        %   Author: Gianluca Iori (gianluca.iori@charite.de)
        %   BSRT - Charite Berlin
        %   Created on:   22/10/2016
        %   Last update:  22/10/2016
        %   ______________________________________________________

        if nargin < 2,   disp('writeMHD: filename not specified: no data will be written..');    return;   end

        %% list sctdata object properties and create parameter structur for MHD output
        p = properties(this);
        for i=1:length(p)
            switch (char(p{i}))
                case 'data'
                    image.data = permute(this.data, [2 1 3]);
                case 'dims'
                    par.ndims = this.dims;
%                     par.ndims = 3;              % required for the use of ITK
                case 'size'
                    par.dimsize = this.size;
                case 'voxelsize'
                    par.elementspacing = this.voxelsize;
                case 'binarydata'
                    if this.binarydata
                        par.binarydata  = 'True';
                    else
                        par.binarydata  = 'False';
                    end
                case 'binarydatabyteordermsb'
                    par.binarydatabyteordermsb = this.binarydatabyteordermsb;
                case 'compresseddata'
                    par.compresseddata = this.compresseddata;
                case 'transformmatrix'
                    if isempty(this.transformmatrix)
                        par.transformmatrix = [1 0 0 0 1 0 0 0 1];
                    else
                        par.transformmatrix = this.transformmatrix;
                    end
                case 'centerofrotation'
                    if isempty(this.centerofrotation)
                        par.centerofrotation = [0 0 0];
                    else
                        par.centerofrotation = this.centerofrotation;
                    end
                case 'offset'
                    if isempty(this.offset)
                        par.offset = [0 0 0];
                    else
                        par.offset = this.offset;
                    end
                case 'anatomicalorientation'
                    par.anatomicalorientation = this.anatomicalorientation;
                case 'elementnumberofchannels'
                    par.elementnumberofchannels = this.elementnumberofchannels;
                case 'elementtype'
                    par.elementtype = this.elementtype;
            end
        end

        %% build argin cell array
        varnames    = fieldnames(par);
        varvalues   = struct2cell(par);
        argin = cell(2*length(varnames),1);
        for i=1:length(varnames)
            argin(2*i-1) = varnames(i);     argin(2*i) = varvalues(i);
        end

        %% write MHD data
        MHDdata.write(filename, image, argin{:});

    end

    function writeSlice(this, filename, mode, slicen, crop)
    %function WRITESLICE(this, filename, mode)
    %   print slice of volume data
    %
    %       mode:   0 only GV data (default)
    %               1 only mask
    %               2 overlap of data and mask (transparent mask; export_fig)
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   29/01/2019
    %   Last update:  29/01/2019
    %
    %   this function is part of the synchro toolbox    
    %   _____________________________________________________
    
    if nargin < 5,          crop = false;    end
    if isempty(crop),   	crop = false;    end
    
    if nargin < 4,          slicen = [];    end
    if isempty(slicen),   	slicen = this.XYslicen;   end

    if nargin < 3,          mode = 0;   end
    if isempty(mode),   	mode = 0;   end

    if nargin < 2,          filename = this.headerfile;   end
    if isempty(filename),   filename = this.headerfile;   end

    if isempty(this.mask),  mode = 0;   end
    
    % check the filename
    [FILEPATH, NAME, EXT] = fileparts(filename);
    NAME = [FILEPATH filesep NAME];
    
    % get the slice and the corresponding mask
    slice   = this.getSlice(slicen, 3);
    if mode>0
        sliceBW = this.getSliceMask(slicen, 3);
    end
    
    % crop both with bounding box
    if crop
        [tmp, level] = threshbone(slice(slice>0));
        [row0, rowd, col0, cold] = bbox(threshbone(slice, 'fixed', level), 10, 0);
        slice = slice(row0:row0+rowd, col0:col0+cold);
        if mode>0
            sliceBW = sliceBW(row0:row0+rowd, col0:col0+cold);
        end
    end
    
    % scale to 8bit unsigned integer
    switch class(slice)
        case {'double' 'single' 'int16' 'uint16'}
            slice = double(slice);
            onepercent = quantile(slice(slice>0), 0.01);
            slice = slice - onepercent;
            slice = uint8((slice./(max(slice(:))))*255);
    end
    
    % print the image
    switch mode
        case 0
            imwrite(slice,[NAME '_s' num2str(slicen) '.tiff'],'Compression','none');
        case 1
            imwrite(sliceBW,[NAME '_BWs' num2str(slicen) '.tiff'],'Compression','none');
        case 2
            % with imshowpair
%             figure; imshowpair(slice, sliceBW); axis image;
%             fig = gcf; export_fig([NAME '_s' slicen], '-dpng'); close(fig);
            
            % with transparent mask
            figure('units','normalized','outerposition',[0 0 1 1]);
            green = cat(3, zeros(size(slice)), ones(size(slice)), zeros(size(slice)));
            h = imagesc(slice);
            colormap gray;
            axis image;
            hold on;
            h = imagesc(green);
            hold off;
            set(h, 'AlphaData', 0.3*double(sliceBW));
            axis off;
            fig = gcf; export_fig([NAME '_s' num2str(slicen)], '-tif', '-native');
            close(fig);
    end
    end
    
    function writeMidplanes(this, filename)
    %function WRITEMIDPLANES(this, filename)
    %   print volume data midplanes
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   03/12/2017
    %   Last update:  08/07/2018
    %
    %   this function is part of the synchro toolbox    
    %   _____________________________________________________
    
    if nargin < 2,   filename = this.headerfile;   end

    switch this.dims
        case 3
            [XY,XZ,YZ] = midplanes(this.data);
            
            switch class(XY)
                case {'double' 'single' 'int16' 'uint16'}
                    XY = XY-min(XY(:));
                    XY = uint8((XY./(max(XY(:))))*255);

                    XZ = XZ-min(XZ(:));
                    XZ = uint8((XZ./(max(XZ(:))))*255);

                    YZ = YZ-min(YZ(:));
                    YZ = uint8((YZ./(max(YZ(:))))*255);
            end
            
            if isempty(this.mask)
                % save pngs
                imwrite(XY,[filename(1:end-4) '_XY.png'],'Compression','none');
                imwrite(XZ,[filename(1:end-4) '_XZ.png'],'Compression','none');
                imwrite(YZ,[filename(1:end-4) '_YZ.png'],'Compression','none');
            else
                [XYm,XZm,YZm] = midplanes(this.mask);
                % save pngs
                figure; imshowpair(XY, XYm); axis image;
                fig = gcf; export_fig([filename(1:end-4) '_XY'], '-dpng'); close(fig);
                figure; imshowpair(XZ, XZm); axis image;
                fig = gcf; export_fig([filename(1:end-4) '_XZ'], '-dpng'); close(fig);
                figure; imshowpair(YZ, YZm); axis image;
                fig = gcf; export_fig([filename(1:end-4) '_YZ'], '-dpng'); close(fig);
            end

            case 2
                XY = this.data;

                switch class(XY)
                    case {'double' 'single' 'int16' 'uint16'}
                        XY = XY-min(XY(:));
                        XY = uint8((XY./(max(XY(:))))*255);
                end
                
                if isempty(this.mask)
                    % save pngs
                    imwrite(XY,[filename(1:end-4) '_XY.png'],'Compression','none');
                else
                    XYm = this.mask;
                    % save pngs
                    figure; imshowpair(XY, XYm); axis image;
                    fig = gcf; export_fig([filename(1:end-4) '_XY'], '-dpng'); close(fig);
                end
    end 
    end

    function writeMATData(this, filename)
    %function WRITEMATDATA(this, filename)
    %   writes volume data (only the data!) as .MAT file
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   04/12/2017
    %   Last update:  04/12/2017
    %
    %   this function is part of the synchro toolbox    
    %   _____________________________________________________
    
    if nargin < 2,   filename = this.headerfile;   end

        % save .MAT
        save([filename(1:end-4) '.mat'], 'this.data', '-v7.3');
        fprintf('data field written to file..\n');

    end

    function writeMAT(this, filename)
    %function WRITEMAT(this, filename)
    %   writes sctdata object to .MAT file
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   04/12/2017
    %   Last update:  04/12/2017
    %
    %   this function is part of the synchro toolbox    
    %   _____________________________________________________
    
    if nargin < 2,   filename = this.headerfile;   end

        % save .MAT
        save([filename(1:end-4) '.mat'], 'this', '-v7.3');
        fprintf('sctdata obj written to file..\n');

    end
    
    %% sctdata setters
    function setMu_scaling(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty','finite'},mfilename,'mu_scaling',1);
        this.mu_scaling = prop;
    end

    function setSlope_dens(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty','finite'},mfilename,'slope_dens',1);
        this.slope_dens = prop;
    end

    function setOffset_dens(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty','finite'},mfilename,'offset_dens',1);
        this.offset_dens = prop;
    end
    
    function setSize(this, sizein)
        validateattributes(sizein,{'numeric'},{'nonempty','finite','positive','integer'},mfilename,'size');
        this.size = sizein;
    end

    function setDims(this, dimsin)
        validateattributes(dimsin,{'numeric'},{'nonempty','finite','positive','integer'},mfilename,'dims',3);
        this.dims = dimsin;
    end

    function setVoxelsize(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty','finite','positive'},mfilename,'voxelsize',1);
        this.voxelsize = prop;
    end
    
    function setOffset(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty','finite'},mfilename,'offset',1);
        this.offset = prop;
    end
    
    function setTransformmatrix(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty','finite'},mfilename,'transformmatrix',1);
        this.transformmatrix = prop;
    end
    
    function setData( this, prop )
    % Setter for numeric property data
    % Usage:
    % 	 this.setData(  prop )
    % Parameters:
    % 	this:  @type 
    % 	prop:  @type 
    
        % validateattributes(prop,{'numeric'},{'nonempty','finite'},mfilename,'data');
        validateattributes(prop,{'numeric', 'logical'},{'nonempty'},mfilename,'data');
        
        this.data = prop;
        fprintf('data set..\n');

        this.setDims(ndims(prop));
        datasize = size(prop);
        this.setSize(datasize);
        
        switch this.dims
            case 3
                this.setXYslicen(round(datasize(3)/2));
                this.setXZslicen(round(datasize(1)/2));
                this.setYZslicen(round(datasize(2)/2));
                this.setOrthoslices;

            case 2
                this.setXYslicen(1);
                this.setXYslice;
        end

    end
    
    function setMask( this, prop )
    % Setter for binary mask mask property
    % Usage:
    % 	 this.setMask(  prop )
    % Parameters:
    % 	this:  @type 
    % 	prop:  @type 
    
        validateattributes(prop,{'logical'},{'nonempty','finite'},mfilename,'mask');
        
        if size(prop) ~= size(this.data)
              error('sctdata.setMask: mask must have same size as data!');
        end
        
        this.mask = prop;
        fprintf('binary data set..\n');

        switch this.dims
            case 3
                this.setXYslicemask;
                this.setXZslicemask;
                this.setYZslicemask;

            case 2
                this.setXYslicemask;
        end
        
    end
    
    function setXYslicen(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty', 'finite', 'positive'},mfilename,'XYslicen');
        this.XYslicen = prop;
    end

    function setXZslicen(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty', 'finite', 'positive'},mfilename,'XZslicen');
        this.XZslicen = prop;
    end

    function setYZslicen(this, prop)
        validateattributes(prop,{'numeric'},{'nonempty', 'finite', 'positive'},mfilename,'YZslicen');
        this.YZslicen = prop;
    end

    function setXYslice(this)
        this.XYslice = this.getSlice(this.XYslicen, 3);
        this.setXYslicemask;
    end
    
    function setXYslicemask(this)
        if ~isempty(this.mask)
            this.XYslicemask = this.getSliceMask(this.XYslicen, 3);
        end
    end
    
    function setYZslice(this)
        this.YZslice = this.getSlice(this.YZslicen, 2);
        this.setYZslicemask;
    end
    
    function setYZslicemask(this)
        if ~isempty(this.mask)
            this.YZslicemask = this.getSliceMask(this.YZslicen, 2);
        end
    end
    
    function setXZslice(this)
        this.XZslice = this.getSlice(this.XZslicen, 1);
        this.setXZslicemask;
    end
    
    function setXZslicemask(this)
        if ~isempty(this.mask)
            this.XZslicemask = this.getSliceMask(this.XZslicen, 1);
        end
    end
    
    function setOrthoslices(this)
        if isempty(this.data)       return;     end
        this.setXYslice;
        this.setXZslice;
        this.setYZslice;
    end
    
    function maxxproj = calcAndSetMax_xproj(this)
        maxxproj = this.calcMax_xproj;
        this.max_xproj = maxxproj;
    end
    
    function maxyproj = calcAndSetMax_yproj(this)
        maxyproj = this.calcMax_yproj;
        this.max_yproj = maxyproj;
    end
    
    function maxzproj = calcAndSetMax_zproj(this)
        maxzproj = this.calcMax_zproj;
        this.max_zproj = maxzproj;
    end
    
    function calcAndSetMax_projs(this)
        this.calcAndSetMax_xproj;
        this.calcAndSetMax_yproj;
        this.calcAndSetMax_zproj;
    end

    %% sctdata getters
    function GetSize(this)
        % size on the physical memory of the computer
        props = properties(this);
        totSize = 0;
        for ii=1:length(props)
            currentProperty = getfield(this, char(props(ii)));
            s = whos('currentProperty'); totSize = totSize + s.bytes;
        end
        fprintf(1, '%d bytes\n', totSize);
    end

    function imsize = get.size(this)
        if ~isempty(this.data)
            if ndims(this.data) == 3
                imsize = [size(this.data,2) size(this.data,1) size(this.data,3)];
            elseif ndims(this.data) == 2
                imsize = [size(this.data,2) size(this.data,1)];
            end
        else
            imsize = this.size;
        end
    end
    
    function imdims = get.dims(this)
        if size(this.data, 3) > 1
            imdims = 3;
        elseif size(this.data, 3) == 1
            imdims = 2;
        else
            imdims = [];
        end
    end

    function sliceim = getSlice(this, sln, dir)
    % returns volume slice at position (slice number) "sln" and among direction 'dir' as 2D matrix 
    %  
    %    if dir = 1 (y) -> sliceim is a this.size(1) x this.size(3) matrix 
    %    if dir = 2 (x) -> sliceim is a this.size(2) x this.size(3) matrix 
    %    if dir = 3 (z) -> sliceim is a this.size(1) x this.size(2) matrix 
    %   ______________________________________________________
    %
    %   Author:         Daniel Rohrbach
    
        validateattributes(sln,{'numeric'},{'nonempty','finite','positive'},mfilename,'sln',1);
        validateattributes(dir,{'numeric'},{'nonempty','finite','positive','>',0,'<',4},mfilename,'dir',2);
        
        switch dir
            case 2
                %permute dimension order to get a well image
                if this.dims < 4
                    d = permute( this.getData(), [3, 1, 2] );
                    sliceim = d(:,:,sln);
                else
                    d = this.getData;               
                    d = permute( d(:,:,:,:), [3, 1, 2, 4] );
                    sliceim = d(:,:,sln);
                end
            case 1
                %permute dimension order to get a well image
                if this.dims < 4
                    d = permute( this.getData(), [3, 2, 1] );
                    sliceim = d(:,:,sln);
                else
                    d = this.getData;
                    d = permute( d(:,:,:,:), [3, 2, 1, 4] );
                    sliceim = d(:,:,sln);
                end

            case 3
                % in z direction we do not need any permutation
                d = this.getData();
                sliceim = d(:,:,sln);
        end
    end
    
    function sliceim = getSliceMask(this, sln, dir)
    % returns volume mask slice at position (slice number) "sln" and among direction 'dir' as 2D matrix 
    %  
    %    if dir = 1 (y) -> sliceim is a this.size(1) x this.size(3) matrix 
    %    if dir = 2 (x) -> sliceim is a this.size(2) x this.size(3) matrix 
    %    if dir = 3 (z) -> sliceim is a this.size(1) x this.size(2) matrix 
    %   ______________________________________________________
    %
    %   Author:         Daniel Rohrbach
    
        validateattributes(sln,{'numeric'},{'nonempty','finite','positive'},mfilename,'sln',1);
        validateattributes(dir,{'numeric'},{'nonempty','finite','positive','>',0,'<',4},mfilename,'dir',2);
        
        switch dir
            case 2
                %permute dimension order to get a well image
                if this.dims < 4
                    d = permute( this.getMask(), [3, 1, 2] );
                    sliceim = d(:,:,sln);
                else
                    d = this.getMask;               
                    d = permute( d(:,:,:,:), [3, 1, 2, 4] );
                    sliceim = d(:,:,sln);
                end
            case 1
                %permute dimension order to get a well image
                if this.dims < 4
                    d = permute( this.getMask(), [3, 2, 1] );
                    sliceim = d(:,:,sln);
                else
                    d = this.getMask;
                    d = permute( d(:,:,:,:), [3, 2, 1, 4] );
                    sliceim = d(:,:,sln);
                end

            case 3
                % in z direction we do not need any permutation
                d = this.getMask();
                sliceim = d(:,:,sln);
        end
    end
    
    function data = getData(this)
        data = this.data;
    end

    function data = getMask(this)
        data = this.mask;
    end

    function minimum = getMin_data_value(this)
         minimum = min(min(min(this.data)));
    end
        
    function maximum = getMax_data_value(this)
         maximum = max(max(max(this.data)));
    end
    
    function minxproj = calcMin_xproj(this)
         minxproj = min(this.data, [], 2);
         minxproj = permute(minxproj,[3 2 1]);
    end

    function maxxproj = calcMax_xproj(this)
        maxxproj = max(this.data, [], 2);
        maxxproj = permute(maxxproj,[3 1 2]);
    end
    
    function minyproj = calcMin_yproj(this)
         minyproj = min(this.data, [], 1);
    end

    function maxyproj = calcMax_yproj(this)
        maxyproj = max(this.data, [], 1);
        maxyproj = permute(maxyproj,[3 2 1]);
    end

    function minzproj = calcMin_zproj(this)
         minzproj = min(this.data, [], 3);
    end

    function maxzproj = calcMax_zproj(this)
        maxzproj = max(this.data, [], 3);
    end

    function datatype = get.elementtype(this)
    % data type
         datatype = class(this.data);
    end
    
    %% sctdata class constructor
    function this = sctdata()
       this.data = []; 
       this.dims = 3;
       this.size = [0 0 0];
       this.voxelsize = [1 1 1];
       this.elementtype = '';
       this.rawfile = '';
       this.headerfile = '';
       this.transformmatrix = [];
       this.centerofrotation = [];
       this.offset = [];
       this.dataoffset = [];
       this.intensity = [];
       this.energy = [];
       this.mu_scaling = [];
       this.slope_dens = [];
       this.offset_dens = [];
       this.nr_of_projections = [];
       this.min_data_value = [];
       this.max_data_value = [];
       this.XYslice = [];
       this.YZslice = [];
       this.XZslice = [];
       this.XYslicemask = [];
       this.YZslicemask = [];
       this.XZslicemask = [];
       this.XYslicen = [];
       this.YZslicen = [];
       this.XZslicen = [];
       this.min_xproj = [];
       this.max_xproj = [];
       this.min_yproj = [];
       this.max_yproj = [];
       this.min_zproj = [];
       this.max_zproj = [];
       this.nr_of_bytes = [];
       this.binarydata = true;
       this.mask = [];
       this.thresholdvalue = [];
       this.BMDdata = false;
    end

end

methods (Static)
%     [result, header] = MHDLoad(filename, slicerange)
end

end
