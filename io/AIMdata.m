classdef AIMdata < handle
    %AIMDATA class for Scanco Image AIM data read
    %
    %   The AIM file extension is associated with the high-resolution computed tomography picture image file format for microCT developed by SCANCO Medical AG.
    %   The *.aim files are viewable in the Bio-Formats software.
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   13/04/2017
    %   Last update:  13/02/2018
    %
    %   this class is part of the synchro toolbox
    %   ______________________________________________________
    
    properties
    end
    
    methods
    end

    methods (Static=true)
        
        function [header, data] = load(filename)
            %LOAD loads Scanco AIM CT data
            %   
            %   data2004_R = AIMLoad('goofy.AIM');
            %   ______________________________________________________
            %
            %   Last update:  10/01/2018
            %   ______________________________________________________

            % get AIM filename if no file is given as input
            if nargin == 0
                [filename, pathname] = uigetfile({'*.aim;*.AIM','Scanco AIM data (*.aim,*.AIM)';'*.*', 'All Files (*.*)'},'Load Scanco AIM data');
                if filename==0
                    % user pressed cancel
                    return;
                end
                
                filename = [pathname filename];
                
            end

            header = AIMdata.readheader(filename, false);

            data = zeros(header.x_dim,header.y_dim,header.z_dim,'int16');

            %% Open file and seek header
            fid = fopen(filename, 'rb');
            fseek(fid, header.offset, -1);

            % 	%% Seek initial frames
            % 	fseek(fid, header.x_dim*header.y_dim*(z_min-1)*2, 0);

            %% read data
            fprintf('Loading AIM data...\n');    

            for z = 1:header.z_dim
                data(:,:,z) = fread(fid,[header.x_dim,header.y_dim],'int16');
            end

            fclose(fid);

            data = permute(data, [2 1 3]);
            fprintf(' done!\n');

        end
        
        function [header, fid] = readheader(filename, leaveopen)
            %READHEADER reads Scanco AIM file header
            %   [header] = AIMdata.readheader(filename)
            %   [header, fid] = AIMReadHeader(filename, 1)  reads header and leaves file open
            %   
            %   The current state of this function is a basic implementation.
            %   Only few crucial information are read from the header!
            %   ______________________________________________________
            %
            %   Last update:  13/04/2017
            %   ______________________________________________________
            
            if ~exist('leaveopen','var') || isempty(leaveopen)
                leaveopen = false;
            end

            %% FOPEN
            fid = fopen(filename,'r');                                          % open file

            if(fid == -1)                                                       % check if opening was successfull
                error(sprintf('Cannot open file %s\n', filename));
            end

            %% READ HEADER
            % first bytes
            offset = fread(fid,3,'uint');
            header.offset = sum(offset);
            % check = fread(fid,4,'char');                                        % char check[16]
            % header.data_type_id = fread(fid,1,'uint');                      % int data_type
            header.type = 'int16';                                          % standard for ISQ
            % header.nr_of_bytes = fread(fid,1,'uint');                       % int nr_of_bytes; /* either one of them */

            % dimensions
            fseek(fid,56,-1);                                                   % skip 56 bytes to x_dim
            header.x_dim = fread(fid,1,'uint');                             % read x_dim [pixels]
            header.y_dim = fread(fid,1,'uint');                             % read y_dim [pixels]
            header.z_dim = fread(fid,1,'uint');                             % read z_dim [slice num]
            % header.x_dim_um = fread(fid,1,'uint');                          % read x_dim_um [um]
            % header.y_dim_um = fread(fid,1,'uint');                          % read y_dim_um [um]
            % header.z_dim_um = fread(fid,1,'uint');                          % read z_dim_um [um]
            % header.slice_thickness_um = fread(fid,1,'uint');                % slice_thickness_um
            % header.slice_increment_um = fread(fid,1,'uint');                % slice_increment_um
            % header.slice_1_pos_um = fread(fid,1,'int');                    % slice_1_pos_um
            % header.min_data_value = fread(fid,1,'int');                    % min_data_value
            % header.max_data_value = fread(fid,1,'int');                    % max_data_value
            % header.mu_scaling = fread(fid,1,'int');                        % mu_scaling (p(x,y,z)/mu_scaling = value [1/cm])
            % header.nr_of_samples = fread(fid,1,'uint');                     % nr_of_samples
            % header.nr_of_projections = fread(fid,1,'uint');                 % nr_of_projections
            % header.scandist_um = fread(fid,1,'uint');                       % scandist_um
            % header.scanner_type = fread(fid,1,'int');                      % scanner_type
            % header.sampletime_us = fread(fid,1,'uint');                     % sampletime_us
            % header.index_measurement = fread(fid,1,'int');                 % index_measurement
            % header.site = fread(fid,1,'int');                              % site /* Coded value */
            % header.reference_line_um = fread(fid,1,'int');                 % reference_line_um
            % header.recon_alg = fread(fid,1,'int');                         % recon_alg /* Coded value */
            % samplename = fread(fid,40,'*char');                                 % read sample name
            % header.samplename = samplename';
            % header.energy = fread(fid,1,'uint');                            % energy [V]
            % header.intensity = fread(fid,1,'uint');                         % intensity [uA]

            % data offset
            % fseek(fid,508,-1);                                                  % fseek to last 4 header bytes
            % offset = fread(fid,1,'uint');                                       % data_offset
            % header.offset = offset * 512 + 512;                             % [bits]
            % header.offset = 3196;
            % header.offset = 3197;                                              % [bytes]

            if leaveopen == false
                fclose(fid);
            end

            end
    end
    
end

