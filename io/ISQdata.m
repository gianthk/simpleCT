classdef ISQdata < handle
    %ISQDATA class for Scanco Medical ISQ data input
    %   Detailed explanation goes here
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Last update:  01/02/2018
    %
    %   this class is part of the synchro toolbox    
    %   ______________________________________________________

    properties
    end
    
    methods
    end
    
    methods (Static=true)
        
        function [data, header, filename] = load(filename, x_min, y_min, z_min, x_size, y_size, z_size, interactive)
            %Load ISQ file
            %   ______________________________________________________
            %
            %   Author:     Vantte Kilappa
            %   Maintainer: Gianluca Iori (gianluca.iori@charite.de)
            %   BSRT - Charite Berlin
            %   Created on:   --/--/----
            %   Last update:  09/01/2018
            %   ______________________________________________________

            if nargin < 8
                interactive = true;         % flag for interactive select of load ranges
            else
              	validateattributes(interactive,   {'logical','numeric'}, {'binary'}, mfilename,'interactive',3);
            end
            
            % get ISQ filename if no file is given as input
            if nargin == 0
                [filename, pathname] = uigetfile({'*.isq;*.ISQ','Scanco ISQ data (*.isq,*.ISQ)';'*.*', 'All Files (*.*)'},'Load Scanco ISQ data');
                if filename==0
                    % user pressed cancel
                    return;
                end
                
                filename = [pathname filename];
            end
            
            % if dims are given as input check input validity
            if nargin == 7
                validateattributes(x_min,   {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'x_min',3);
                validateattributes(y_min,   {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'y_min',3);
                validateattributes(z_min,   {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'z_min',3);
                validateattributes(x_size,  {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'x_size',3);
                validateattributes(y_size,  {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'y_size',3);
                validateattributes(z_size,  {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'z_size',3);
                interactive = false;
            end
            
            %  throw warning if the call was not correct
            if nargin>1 && nargin<7
                warning('ISQdata.load: Not enough inputs.. Example call: [data, header, filename] = ISQdata.load(filename, x_min, y_min, z_min, x_size, y_size, z_size)');
            end
            
            % fetch ISQ header info
            header = ISQdata.readheader(filename);
            
            % interactive select load range if not given as input
            if interactive
                
                % dialog box to select load range 
                fprintf('select portion of ISQ data to load..\n');
                prompt={'Xmin [voxels]','Ymin [voxels]','Zmin [voxels]','Xsize [voxels]','Ysize [voxels]','Zsize [voxels]'};
                dlg_title = 'ISQ data load range';
                defaultans = {'1' '1' '1' num2str(header.x_dim) num2str(header.y_dim) num2str(header.z_dim)};
                num_lines = [1 14; 1 14; 1 14; 1 14; 1 14; 1 14];
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                
                x_min = str2num(answer{1,:});
                y_min = str2num(answer{2,:});
                z_min = str2num(answer{3,:});
                x_size = str2num(answer{4,:});
                y_size = str2num(answer{5,:});
                z_size = str2num(answer{6,:});

                validateattributes(x_min,   {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'x_min',3);
                validateattributes(y_min,   {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'y_min',3);
                validateattributes(z_min,   {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'z_min',3);
                validateattributes(x_size,  {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'x_size',3);
                validateattributes(y_size,  {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'y_size',3);
                validateattributes(z_size,  {'numeric'}, {'nonnegative','integer','nonempty','finite'}, mfilename,'z_size',3);
                
            end
            
            % read data
            data = ISQdata.readdata(filename, x_min, y_min, z_min, x_size, y_size, z_size);
            
        end
        
        function [headerinfo, fid] = readheader(filename, leaveopen)
            %ISQREADHEADER reads Scanco ISQ file header
            %   [headerinfo] = ISQdata.readheader(filename,leaveopen)
            %   [headerinfo,fid] = ISQdata.readheader(filename,leaveopen)  reads header and leaves file open
            %   ______________________________________________________
            %
            %   Author: Gianluca Iori (gianluca.iori@charite.de)
            %   BSRT - Charite Berlin
            %   Created on:   16/11/2015
            %   Last update:  09/01/2018
            %   ______________________________________________________
            
            if ~exist('leaveopen','var') || isempty(leaveopen)
                leaveopen = false;
            end

            % FOPEN
            fid = fopen(filename,'r');                                          % open file

            if(fid == -1)                                                       % check if opening was successfull
                error(sprintf('Cannot open file %s\n',filename));
            end

            % READ HEADER
            % first bytes
            check = fread(fid,4,'char');                                        % char check[16]
            headerinfo.data_type_id = fread(fid,1,'uint');                      % int data_type
            headerinfo.type = 'int16';                                          % standard for ISQ
            headerinfo.nr_of_bytes = fread(fid,1,'uint');                       % int nr_of_bytes; /* either one of them */

            % dimensions
            fseek(fid,44,-1);                                                   % skip 44 bytes to x_dim
            headerinfo.x_dim = fread(fid,1,'uint');                             % read x_dim [pixels]
            headerinfo.y_dim = fread(fid,1,'uint');                             % read y_dim [pixels]
            headerinfo.z_dim = fread(fid,1,'uint');                             % read z_dim [slice num]
            headerinfo.x_dim_um = fread(fid,1,'uint');                          % read x_dim_um [um]
            headerinfo.y_dim_um = fread(fid,1,'uint');                          % read y_dim_um [um]
            headerinfo.z_dim_um = fread(fid,1,'uint');                          % read z_dim_um [um]
            headerinfo.slice_thickness_um = fread(fid,1,'uint');                % slice_thickness_um
            headerinfo.slice_increment_um = fread(fid,1,'uint');                % slice_increment_um
            headerinfo.slice_1_pos_um = fread(fid,1,'int');                    % slice_1_pos_um
            headerinfo.min_data_value = fread(fid,1,'int');                    % min_data_value
            headerinfo.max_data_value = fread(fid,1,'int');                    % max_data_value
            headerinfo.mu_scaling = fread(fid,1,'int');                        % mu_scaling (p(x,y,z)/mu_scaling = value [1/cm])
            headerinfo.nr_of_samples = fread(fid,1,'uint');                     % nr_of_samples
            headerinfo.nr_of_projections = fread(fid,1,'uint');                 % nr_of_projections
            headerinfo.scandist_um = fread(fid,1,'uint');                       % scandist_um
            headerinfo.scanner_type = fread(fid,1,'int');                      % scanner_type
            headerinfo.sampletime_us = fread(fid,1,'uint');                     % sampletime_us
            headerinfo.index_measurement = fread(fid,1,'int');                 % index_measurement
            headerinfo.site = fread(fid,1,'int');                              % site /* Coded value */
            headerinfo.reference_line_um = fread(fid,1,'int');                 % reference_line_um
            headerinfo.recon_alg = fread(fid,1,'int');                         % recon_alg /* Coded value */
            samplename = fread(fid,40,'*char');                                 % read sample name
            headerinfo.samplename = samplename';
            headerinfo.energy = fread(fid,1,'uint');                            % energy [V]
            headerinfo.intensity = fread(fid,1,'uint');                         % intensity [uA]

            % data offset
            fseek(fid,508,-1);                                                  % fseek to last 4 header bytes
            offset = fread(fid,1,'uint');                                       % data_offset
            headerinfo.offset = offset * 512 + 512;                             % [bits]

            if leaveopen == false
                fclose(fid);
            end
            
        end
        
        function data = readdata(filename, x_min, y_min, z_min, x_size, y_size, z_size)
            %LOADDATA read portion of Scanco ISQ data
            %   
            %   data2004_R = ISQdata.readdata('S:\#Charite-Central\BCRT\AG-raum-qbam-archiv-read\2015.003.qbam.TaCoSound\data\XtremeCT-II\femur\2004_R\C0001898.ISQ',1,1,400,4600,4600,1);
            %   ______________________________________________________
            %
            %   Author:     Vantte Kilappa
            %   Maintainer: Gianluca Iori (gianthk.iori@gmail.com)
            %   BSRT - Charite Berlin
            %   Created on:     --/--/----
            %   Last update:    09/01/2018
            %   ______________________________________________________
            
            if ~exist(filename, 'file') 
                error('Input file missing!');
            end

            headerinfo = ISQdata.readheader(filename, false);

            data = zeros(x_size,y_size,z_size,'int16');

            % Open file and seek header
            fid = fopen(filename, 'rb');
            fseek(fid, headerinfo.offset, -1);

            % Seek initial frames
            fseek(fid, headerinfo.x_dim*headerinfo.y_dim*(z_min-1)*2, 0);

            fprintf('Reading ISQ data...\n');    
            % For loop for reading data
            for kk=1:z_size
        %         fprintf('kk=%i\n',kk);
                % Seek everything before the first line we want to read
                fseek(fid, headerinfo.x_dim*(y_min-1)*2, 0);
                for jj=1:y_size
        %             fprintf('jj=%i\n',jj);
                    % Seek everything before the first column we want to read
                    fseek(fid, (x_min-1)*2, 0);
                    % Read line
                    data(:, jj, kk) = fread(fid, x_size, 'int16');
                    % Seek everything after the last column we want to read
                    fseek(fid, (headerinfo.x_dim-(x_min+x_size-1))*2, 0);
                end
                % Seek everything after the last line we want to read
                fseek(fid, headerinfo.x_dim*(headerinfo.y_dim-(y_min+y_size-1))*2, 0);
            end
            fclose(fid);
            fprintf(' done!\n');
            
            % compensate for matlab row-column convention
            data = permute(data, [2 1 3]);

        end
        
        function data = readslices(filename, offset, rows, cols, zmin, zmax)
            %READSLICES fast read for (slices) portion of ISQ data
            %   ______________________________________________________
            %
            %   Author:         Vantte Kilappa
            %   Maintainer:     Gianluca Iori (gianthk.iori@gmail.com)
            %   BSRT - Charite Berlin
            %   Created on:   --/--/2015
            %   Last update:  10/01/2018
            %   ______________________________________________________
            
            fid = fopen(filename, 'rb');
            fseek(fid, offset, -1);
            fseek(fid, rows*cols*(zmin-1)*2, 0);
            zdim = zmax-zmin+1;
            data = fread(fid, rows*cols*zdim, 'int16');
            data = reshape(data,rows,cols,zdim);
            data = permute(data, [2 1 3]);              % compensate for matlab row-column convention
            fclose(fid);
        end
        
        function data = resample(filename, resamplefactor, method)
            %RESAMPLE loads and resamples large ISQ data to smaller size
            %   data = ISQDATA.RESAMPLE(pippo.ISQ, RESAMPLEFACTOR, METHOD) resamples input ISQ data with given RESAMPLEFACTOR and METHOD.
            %   Possible values for METHOD are:
            %
            %     'box'         (Default) Averaging box kernel.
            %     'imresize'    Uses MAtlab's IMRESIZE function (not implemented yet).
            %
            %   The 'imresize' METHOD, allows you to specify additional parameters for the control of the actual imresize algorithm.
            %   ______________________________________________________
            %
            %   Authors:        Gianluca Iori, Vantte Kilappa
            %   Maintainer:     Gianluca Iori (gianthk.iori@gmail.com)
            %   BSRT - Charite Berlin
            %   Created on:   --/--/2015
            %   Last update:  10/01/2018
            %   ______________________________________________________
            
            if nargin<3,    method = 'box'; end
            
            method = validatestring(method,{'box','imresize'},mfilename,'METHOD',2);

            switch method
                case 'box'
                    % fetch header info
                    [headerinfo] = ISQdata.readheader(filename);

                    % image matrix size
                    xdim = headerinfo.x_dim;
                    ydim = headerinfo.y_dim;
                    zdim = headerinfo.z_dim;

                    % image matrix size after resample
                    x_last = floor(xdim/resamplefactor);
                    y_last = floor(ydim/resamplefactor);
                    z_last = floor(zdim/resamplefactor);
                    % 	z_last = 521;

                    % initialize resampled data array
                    data = zeros(y_last, x_last, z_last, 'int16');

                    fprintf('ISQdata.resample:\tstart resample file: %s\n', filename);
                    fprintf('\tImage size (original):\t%i x %i x%i\n', xdim, ydim, zdim);
                    fprintf('\tImage size (resample):\t%i x %i x%i\n', x_last, y_last, z_last);

                    tic
                    fprintf('slice:        -/-');                   % complete status output
                    for kk=1:z_last
                        % define read limits
                        z_start = (kk-1)*resamplefactor+1;          % first read slice
                        z_end = kk*resamplefactor;                  % last read slice

                        % read data
                        data_to_resample = ISQdata.readslices(filename, headerinfo.offset, ydim, xdim, z_start, z_end);

                        fprintf(1,[repmat('\b',1,numel(num2str(kk))+numel(num2str(z_last))) '\b%i/%i'], kk, z_last);        % update status
                        for ii=1:y_last
                            y_start = (ii-1)*resamplefactor+1;
                            y_end = ii*resamplefactor;
                            for jj=1:x_last
                                x_start = (jj-1)*resamplefactor+1;
                                x_end = jj*resamplefactor;

                                % resample the data
                                data(ii,jj,kk) = mean(mean(mean(data_to_resample(y_start:y_end,x_start:x_end,:))));
                                % (tested) this is on average faster than tmp = data_to_resample(y_start:y_end,x_start:x_end,:); data(ii,jj,kk)=mean(tmp(:));
                            end
                        end
                    end
                    fprintf('\n')
                    fprintf(' done!\t');
                    toc
                    fprintf('\n')
            
                case 'imresize'
                    error('ISQdata.resample: imresize method not implemented yet!');
            end
        end
       
    end
    
end
