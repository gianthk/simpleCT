classdef DICOMdata < handle
    %basic class for DICOM data input
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   --/--/2015
    %   Last update:  20/05/2018
    %
    %   this class is part of the synchro toolbox
    %   ______________________________________________________
    
    properties
    end
    
    methods
    end

    methods (Static=true)
        
        function [dicomdata, info, FILENAME] = load(filename)
            %Load DICOM file(s)
            %   filename can be cell array of strings containing single slice filenames
            %   example call:   [result, header] = DICOMdata.load('./pippo.DCM');
            %   ______________________________________________________
            
            currentdir = pwd;

            % get DICOM filenames
            if nargin == 0
                [FILENAME, PATHNAME] = uigetfile({'*.DCM;*.dcm','DICOM Files (*.dcm,*.DCM)';'*.*', 'All Files (*.*)'},'Load DICOM Data Files','MultiSelect','on');
                if isequal(FILENAME,0) || isequal(PATHNAME,0)
                    % user pressed cancel
                    return;
                end
                
            else
                if iscell(filename)
                    % cell array of slicename list given
                    [PATHNAME, NAME, EXT] = fileparts(char(filename{1}));
                    for i=1:length(filename)
                        [PATHSTR, NAME, EXT] = fileparts(char(filename{i}));
                        FILENAME{i} = [NAME EXT];
                    end
                else
                    % one single filename given
                    [PATHNAME, NAME, EXT] = fileparts(char(filename));
                    FILENAME{1} = [NAME EXT];
                end
                
                
            end

            % create cell array 
            if isstr(FILENAME)      FILENAMEtmp = FILENAME;     clear FILENAME;     FILENAME{1} = FILENAMEtmp;      end

            % load DICOM data
            fprintf('Loading DICOM data...');
            if iscell(FILENAME)    
                if ischar(PATHNAME) && ~isempty(PATHNAME)
                    % cd file location
                    cd(PATHNAME);
                end

                % get DICOM info
                S = squeeze(dicomread(FILENAME{1}));        % load one image to get data type
                info = dicominfo(FILENAME{1});
                info.Slices = numel(FILENAME);

                % initialize dicomdata
                tp = class(S);
                dicomdata = zeros([info.Rows info.Columns info.Slices], tp);
                % data_tmp = dicomdata;
                % data_tmp = zeros([info.Columns info.Rows info.Slices], tp);
                
                position = zeros(info.Slices,2);

                % load each slice and its position in the stack
                for i=1:info.Slices
                    % data_tmp(:,:,i) = squeeze(dicomread(FILENAME{i}));
                    dicomdata(:,:,i) = squeeze(dicomread(FILENAME{i}));
                    
                    info_tmp = dicominfo(FILENAME{i});
                    if isfield(info_tmp, 'ImagePositionPatient')
                        position(i,:) = [info_tmp.ImagePositionPatient(3) i];
                    else
                        position(i,:) = [i i];
                    end
                end

                % resort the slices according to the image position
                position2 = sortrows(position,1);
                if isequal(position2, position)
                    % directly permute if the position is not messed up
                    % dicomdata = permute(data_tmp, [2 1 3]);
                    dicomdata = permute(dicomdata, [2 1 3]);
                else
                    for i=1:info.Slices
                        % go through position array if slice were load not in the right order
                        % dicomdata(:,:,i) = data_tmp(:,:,position(i,2));
                        % dicomdata(:,:,i) = data_tmp(:,:,position(i,2))';      % compensate for matlab row-column convention
                        dicomdata(:,:,i) = dicomdata(:,:,position(i,2))';      % compensate for matlab row-column convention
                    end
                end
                % clear data_tmp

                cd(currentdir);
                fprintf(' done!\n');
            else
                error('invalid filename');
            end

        end
        
    end
    
end
