classdef imagestackdata < handle
    %basic class for reading stack of images
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   --/--/2015
    %   Last update:  17/10/2017
    %
    %   see also IMREAD, READNRECON_LOG
    %
    %   this class is part of the synchro toolbox
    %   ______________________________________________________
    
    properties
    end
    
    methods
    end

    methods (Static=true)
        
        function [data, info, FILENAME] = load(filename)
            %Load Image Sequence
            %   filename can be cell array of strings containing single slice filenames
            %   example call:   [data, info] = imagestackdata.load('./pippo.png');
            %   ______________________________________________________
            
            currentdir = pwd;

            % get filenames
            if nargin == 0
                [FILENAME, PATHNAME] = uigetfile({'*.jpg;*.gif;*.png;*.bmp;*.tif','Image Files';'*.*','All Files' },'Load Stack of Images','MultiSelect','on');
                if FILENAME==0
                    % user pressed cancel
                    return;
                end
                
            else
                [PATHNAME, NAME, EXT] = fileparts(char(filename{1}));
                for i=1:length(filename)
                    [PATHSTR, NAME, EXT] = fileparts(char(filename{i}));
                    FILENAME{i} = [NAME EXT];
                end
            end

            % create cell array 
            if isstr(FILENAME)      FILENAMEtmp = FILENAME;     clear FILENAME;     FILENAME{1} = FILENAMEtmp;      end

            % load image data
            fprintf('Loading Image Sequence...');
            if iscell(FILENAME)    
                if ischar(PATHNAME) && ~isempty(PATHNAME)
                    % cd file location
                    cd(PATHNAME);
                end

                % load first image of stack
                S = squeeze(imread(FILENAME{1}));
                tp = class(S);              % get data type

                % try to fetch information from NRecon log file
                % we usually load sequencies of microCT data reconstructed with NRecon
                [info] = imagestackdata.readNRecon_log(FILENAME{1});
                info.elementtype = tp;

                % data size
                info.Rows = size(S,1);
                info.Columns = size(S,2);
                info.Slices = numel(FILENAME);

                % pre-allocate data
                data = zeros([info.Rows info.Columns info.Slices], info.elementtype);
                data_tmp = zeros([info.Rows info.Columns info.Slices], info.elementtype);

                % load slices
                position = zeros(info.Slices, 2);                                       % initialize slice position array
                if info.Slices > 1
                    for i=1:info.Slices
                        % fprintf('.');
                        if ndims(squeeze(imread(FILENAME{i}))) > 2
                            data_tmp(:,:,i) = rgb2gray(squeeze(imread(FILENAME{i})));	% load image slice and convert to grayscale if RGB
                        else
                            data_tmp(:,:,i) = squeeze(imread(FILENAME{i}));             % load image slice
                        end
                        position(i,:) = [str2num(FILENAME{i}(end-7:end-4)) i];          % NRecon builds filename leaving slicenumber in last 4 slots
                    end

                    % resort the slices according to the image position
                    position = sortrows(position,1);
                    for i=1:info.Slices
                       data(:,:,i) = data_tmp(:,:,position(i,2)); 
                    end
                else
                    if ndims(squeeze(imread(FILENAME{1}))) > 2
                        data = rgb2gray(squeeze(imread(FILENAME{1})));              % load image slice and convert to grayscale if RGB
                    else
                        data = squeeze(imread(FILENAME{1}));                        % load image slice
                    end
                end
                clear data_tmp

                cd(currentdir);
            else
                error('invalid filename');
            end
            fprintf(' done!\n');
        end
        
        function [info] = readNRecon_log(filename)
            %READNRECON_LOG reads information from NRecon CT data reconctruction log file
            %   [info] = readNRecon_log(filename)
            %   Reads info from reconstruction log file of NRecon CT image data.
            %   Filename of one reconstruction image is taken as input.
            %   ______________________________________________________
            %
            %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
            %   BSRT - Charite Berlin
            %   Created on:   09/12/2015
            %   Last update:  09/12/2015
            %
            %   See also TEXTREAD, MATLAB.LANG.MAKEVALIDNAME.
            %   ______________________________________________________
            
            logfilename = [filename(1:findstr(filename,'rec')+2) '.log'];                       % log file name
            % dims = [0 0 0];

            if exist(logfilename,'file') == 2
                % Read relevant details from log (to be completed)
                [meta_var, meta_val] = textread(logfilename, '%s%s', 'delimiter', '=');

                meta_var = regexprep(meta_var,'\W','');                                         % remove invalid characters from cell fields
            %     meta_var = matlab.lang.makeValidName(meta_var,'ReplacementStyle','delete');     % remove invalid characters from cell fields (Introduced in MATLAB 2014b)
                j=1;
                while j<=length(meta_var)
                    if isvarname(char(meta_var(j)))==0,
                        meta_var(j) = [];   meta_val(j) = [];                                   % delete invalid field names
                    end
                    j=j+1;
                end
                [meta_var,ia,ic] = unique(meta_var);                                            % remove duplicates of fields

                % create info structure with fields and values from log file
                info = cell2struct(meta_val(ia),meta_var,1);

                fields = fieldnames(info);
                for i=1:numel(fields)
                    if all(ismember(info.(fields{i}), '0123456789+-.eEdD'))                     % if field can be converted from string to num
                        info.(fields{i}) = str2num(info.(fields{i}));                           % convert field to numeric
                    end
                end
            else
                warning('NRecon LOG file not found. INFO struct will be empty..');
                info = struct;
            end
        end
    
    end

end
