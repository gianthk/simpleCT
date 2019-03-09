classdef MHDdata < handle
    %MHDDATA class for Medical Image MHD data input/output
    %   ______________________________________________________
    %   Example Header file (.mhd) for a meta image (ITK):
    %             ObjectType = Image
    %             NDims = 3
    %             BinaryData = True
    %             BinaryDataByteOrderMSB = False
    %             CompressedData = False
    %             TransformMatrix = 1 0 0 0 1 0 0 0 1
    %             Offset = 0 0 0
    %             CenterOfRotation = 0 0 0
    %             AnatomicalOrientation = LPS
    %             ElementSpacing = 0.1 0.1 0.1
    %             DimSize = 50 40 30
    %             ElementType = MET_UCHAR
    %             ElementDataFile = coordinates_50_40_30.raw
    %
    %   visit also:
    %       https://de.mathworks.com/matlabcentral/fileexchange/29344-read-medical-data-3d
    %       https://itk.org/Wiki/ITK/MetaIO/Documentation
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:   --/02/2016
    %   Last update:  06/06/2018
    %
    %   this class is part of the synchro toolbox
    %   ______________________________________________________
    
    properties
    end
    
    methods
    end

    methods (Static=true)
        
        function [result, header, filename_mhd] = load(filename, slicerange)
            %Load MHD file
            %   filename must be passed without suffix
            %   ______________________________________________________
            %
            %   Author:     Vantte Kilappa
            %   Maintainer: Gianluca Iori (gianluca.iori@charite.de)
            %   BSRT - Charite Berlin
            %   Created on:   --/--/----
            %   Last update:  10/02/2018
            %   ______________________________________________________

            % get MHD filename if no file is given as input
            if nargin == 0      filename = '';      end
            
            if isempty(filename)
                [filename, pathname] = uigetfile({'*.mhd;*.MHD','MetaImage data (*.mhd,*.MHD)';'*.*', 'All Files (*.*)'},'Load MetaImage data');
                if filename==0
                    % user pressed cancel
                    return;
                end
                
                filename_mhd = [pathname filename];
                filename = [pathname filename(1:end-4)];
            else
                [pathname, filename, ext] = fileparts(filename);
                if ~isempty(pathname)
                    filename = [pathname filesep filename];
                end
            end

            % raw filename
            filename_raw = [filename '.raw'];

            % load data
            fprintf('Loading MHD data...');

            % Fetch info from header
            [header, dims, datatype] = MHDdata.readheader(filename);

            % We load exactly the desired amount of slices and only the desired
            % slices. If slicerange is omitted, we load everything.
            if(~exist('slicerange','var'))
                slicerange = 1:dims(3);
            end
            
            s = dir(filename_raw);          % fetch file info
            
            if ~isempty(s)                  % strange behaviour with '.raw' files not seen by command dir
                if s.bytes == 0
                    % empty file
                    warning('Raw file empty. No data read!');
                    result = [];
                    return;
                end
            end
            
            % open the file and read it
            FILE_IN = fopen(filename_raw, 'rb');

            switch (datatype)
                case 'MET_FLOAT'
                    readmode = 'float32';
                    createmode = 'single';
                    seek_multi = 4;
                case 'MET_DOUBLE'
                    readmode = 'float64';
                    createmode = 'double';
                    seek_multi = 8;
                case 'MET_SHORT'
                    readmode = 'int16';
                    createmode = 'int16';
                    seek_multi = 2;
                case 'MET_USHORT'
                    readmode = 'uint16';
                    createmode = 'uint16';
                    seek_multi = 2;
                case 'MET_UCHAR'
                    readmode = 'uint8';
                    createmode = 'uint8';
                    seek_multi = 1;
                otherwise
                    disp('Unsupported datatype!')
            end

            result = zeros(dims(2),dims(1),length(slicerange), createmode);

            for ii=1:length(slicerange)
                fseek(FILE_IN, (slicerange(ii)-1)*dims(1)*dims(2)*seek_multi, 'bof'); % multiplier 2 due to 16-bit images
                tmp = fread(FILE_IN, [dims(1), dims(2)], readmode);
                result(:,:,ii) = tmp';   % compensate for matlab row-column convention
            end

            fclose(FILE_IN);
            
            fprintf(' done!\n');

        end
        
        function [header, dims, datatype] = readheader(filename, stackheight)
            %read MHD file header
            %   Reads the header of an MHD file. If the optional argument is set,
            %   this returns a dims -line with the selected stack height.
            %   Filename must be passed without the .mhd suffix
            %   If stackheight is passed, the function returns a custom stack height
            %   instead of the real one.
            %
            %   Example Header file (.mhd) for a meta image (ITK):
            %             ObjectType = Image
            %             NDims = 3
            %             BinaryData = True
            %             BinaryDataByteOrderMSB = False
            %             CompressedData = False
            %             TransformMatrix = 1 0 0 0 1 0 0 0 1
            %             Offset = 0 0 0
            %             CenterOfRotation = 0 0 0
            %             AnatomicalOrientation = LPS
            %             ElementSpacing = 0.1 0.1 0.1
            %             DimSize = 50 40 30
            %             ElementType = MET_UCHAR
            %             ElementDataFile = coordinates_50_40_30.raw
            %
            %   visit also:
            %       https://de.mathworks.com/matlabcentral/fileexchange/29344-read-medical-data-3d
            %       https://itk.org/Wiki/ITK/MetaIO/Documentation
            %   ______________________________________________________

            if ~strcmp(filename(end-3:end), '.mhd')
                filename_header = [filename '.mhd'];
            else
                filename_header = filename;
            end

            if(~exist(filename_header, 'file'))
                disp('Header file missing!');
                return;
            end

            dims = [1 1 1];
            ndims = 3;

            % Read relevant details about the raw
            %FILE_IN = fopen(filename_header, 'rt');
            [meta_var, meta_val] = textread(filename_header, '%s%s', 'delimiter', '=');
            header = '';
            for ii=1:length(meta_var)
                header = [header char(meta_var(ii))];
                switch(strtrim(char(meta_var(ii))))
                    case 'NDims'
                        ndims = strread(strtrim(char(meta_val(ii))), '%d');
                        header = [header sprintf('= %s\n', char(meta_val(ii)))];
                    case 'ElementType'
                        datatype = strtrim(char(meta_val(ii)));
                        header = [header sprintf('= %s\n', char(meta_val(ii)))];
                    case 'DimSize'
                        switch ndims
                            case 2
                                [dims(1), dims(2)] = strread(strtrim(char(meta_val(ii))), '%d%d');
                                header = [header sprintf('= %i %i\n', dims(1), dims(2))];
                            case 3
                                [dims(1), dims(2), dims(3)] = strread(strtrim(char(meta_val(ii))), '%d%d%d');
                                if (~exist('stackheight', 'var'))
                                    stackheight = dims(3);
                                end
                                header = [header sprintf('= %i %i %i\n', dims(1), dims(2), stackheight)];
                        end
                    otherwise
                        header = [header sprintf('= %s\n', char(meta_val(ii)))];
                end
            end
        end
        
        function f = write(filename, image, varargin)
            %write ITK meta image data (.mhd)
            %   MHDWrite(filename,image)    writes image to filename
            %   MHDWrite(filename,image,param,value,...)
            %       % params and default value:
            %       'NDims'                     =   3
            %       'BinaryData'                =   true
            %       'BinaryDataByteOrderMSB'    =   false
            %       'CompressedData'            =   false
            %       'TransformMatrix'           =   1 0 0 0 1 0 0 0 1 (coherent with NDims)
            %       'CenterOfRotation'          =   0 0 0 (coherent with NDims)
            %       'AnatomicalOrientation'     =   RAI
            %       'ElementNumberOfChannels'   =   1
            %       'ElementType'               =   MET_FLOAT
            %
            %   Example
            %   -------
            %   save sample 2D image
            %   A = rgb2gray(importdata('ngc6543a.jpg'));
            %   image.size = [size(A,2) size(A,1) 1];
            %   image.orientation = [1 0 0; 0 1 0; 0 0 1];
            %   image.origin = [0 0 0];
            %   image.spacing = [1 1 1];
            %   image.data = A';
            %   write('./pippo.mhd',image);
            %   ______________________________________________________
            %
            %   Author:         Vantte Kilappa
            %   Maintainer:     Gianluca Iori (gianluca.iori@charite.de)
            %   BSRT - Charite Berlin
            %   Created on:   --/02/2016
            %   Last update:  13/03/2018
            %   ______________________________________________________

            if isunix,  filename = strrep(filename,'\','/');    end

            %% element_types struct
            % revise this
            %element_types=struct('uchar','MET_UCHAR','double','MET_DOUBLE','float','MET_FLOAT','logical','MET_CHAR','int8','MET_UCHAR','uint8','MET_UCHAR','int16','MET_SHORT','uint16','MET_USHORT','int32','MET_INT','uint32','MET_UINT');
            element_types = struct('double','MET_DOUBLE','int8','MET_CHAR','uint8','MET_UCHAR','int16','MET_SHORT','uint16','MET_USHORT','int32','MET_INT','uint32','MET_UINT','single','MET_FLOAT','logical','MET_UCHAR');
            element_names = fieldnames(element_types);

            %% Look for default arguments
            if isfield(image,'size'),           NDims                   = numel(image.size);
                                                if numel(image.size) == 2
                                                    DimSize                 = num2str([image.size(:)' 1]);
                                                elseif image.size(3) == 1
                                                    DimSize                 = num2str([image.size(:)' 1]);
                                                else
                                                    DimSize                 = num2str(image.size(:)');
                                                end
                                                
                                                CenterOfRotation        = num2str(zeros(1,NDims));      end

                                                BinaryData              = 'True';
                                                BinaryDataByteOrderMSB  = 'False';
                                                CompressedData          = 'False';
                                                AnatomicalOrientation   = 'RAI';
                                                %ElementType =getfield(element_types,class(image.data));
                                                EType                   = 'double';
                                                ObjectType              = 'image';

            if isfield(image,'origin'),         Offset                  = num2str(image.origin(:)');    end
            if isfield(image,'spacing'),        ElementSpacing          = num2str(image.spacing(:)');   end

            if isfield(image,'orientation'),    TransformMatrix  = num2str(reshape(image.orientation,1,NDims*NDims));   end     %'1 0 0 0 1 0 0 0 1';% (coherent with NDims)

            if (isa(image,'VectorImageType') || isfield(image,'datax'))
                                                ElementNumberOfChannels = 3;
            else                                ElementNumberOfChannels = 1;
            end

            %% Read user specified arguments
            for i=1:size(varargin,2)
                if ischar(varargin{i}),
                    switch (lower(varargin{i}))
                        case 'ndims'
                            NDims = varargin{i+1};
                        case 'binarydata'
                            BinaryData  = varargin{i+1};
                        case 'binarydatabyteordermsb'
                            BinaryDataByteOrderMSB = varargin{i+1};
                        case 'compresseddata'
                            CompressedData = varargin{i+1};
                        case 'transformmatrix'
                            TransformMatrix = varargin{i+1};
                        case 'centerofrotation'
                            CenterOfRotation = varargin{i+1};
                        case 'anatomicalorientation'
                            AnatomicalOrientation = varargin{i+1};
                        case 'elementnumberofchannels'
                            ElementNumberOfChannels = varargin{i+1};
                        case 'elementtype'
                            EType = varargin{i+1};
                        case 'objecttype'
                            ObjectType = varargin{i+1};
                        case 'offset'
                            Offset = varargin{i+1};
                        case 'elementspacing'
                            ElementSpacing = varargin{i+1};
                        case 'dimsize'
                            DimSize = varargin{i+1};
                    end
                end
            end

            % 3rd element of DimSize for 2D images must be 1
            if numel(DimSize) == 2 && NDims == 3
                DimSize(3) = 1;
            end
            
            % if for some reason NDims was set to 2 but data is 3D
            if numel(DimSize) == 3 && NDims == 2
                NDims = 3;
            end
            
            
            if      isfield(element_types,EType) == 0,  EType = element_names(find(strcmp(struct2cell(element_types),EType)));  end
            ElementType =  element_types.(char(EType) );

            % Extract path and filenames---
            path_ = regexp(filename,filesep,'split');
            path = [];
            for i=2:(size(path_,2)-1)
                path = [path filesep path_{i}];
            end

            raw_ = regexp(filename,'\.','split');
            raw = [];
            for i=1:(size(raw_,2)-1)
                raw = [raw  raw_{i} '.'];
            end
            rawfile = [raw 'raw'];
            rawfilename_ = regexp(rawfile,filesep,'split');
            rawfilename = rawfilename_{size(rawfilename_,2)};

            % write *.mhd file
            fid=fopen(filename,'w','native');
            fprintf(fid,'ObjectType = %s\n',ObjectType);
            fprintf(fid,'NDims = %d\n',NDims);
            fprintf(fid,'BinaryData = %s\n',BinaryData);
            fprintf(fid,'BinaryDataByteOrderMSB = %s\n',BinaryDataByteOrderMSB);
            fprintf(fid,'CompressedData = %s\n',CompressedData);
            fprintf(fid,'TransformMatrix = %s\n',num2str(TransformMatrix));
            fprintf(fid,'CenterOfRotation = %s\n',num2str(CenterOfRotation));
            fprintf(fid,'AnatomicalOrientation = %s\n',AnatomicalOrientation);
            fprintf(fid,'Offset = %s\n', num2str(Offset));
            fprintf(fid,'ElementSpacing = %s\n', num2str(ElementSpacing));
            fprintf(fid,'DimSize = %s\n', num2str(DimSize));
            fprintf(fid,'ElementNumberOfChannels = %d\n',ElementNumberOfChannels);
            fprintf(fid,'ElementType = %s\n',ElementType);
            fprintf(fid,'ElementDataFile = %s\n',rawfilename);

            fclose(fid);

            % write *.raw file

            %dataToWrite=image.data;

            %if (islogical(dataToWrite))
            %     dataToWrite = cast(dataToWrite,'uint8');
            % else

            fid=fopen(rawfile,'w','native');
            if (ElementNumberOfChannels ==1)
                %fwrite(fid,image.data,class(dataToWrite) );
                dataAll = cast(image.data,char(EType));    
            elseif (ElementNumberOfChannels==2)
                %do nothing
            elseif (ElementNumberOfChannels==3)
                %write point per point
                dataAll = cast(permute(cat(NDims+1,image.datax,  image.datay, image.dataz),...
                    [NDims+1 1:NDims]),char(EType));
            elseif (ElementNumberOfChannels==4)
                %TODO
            end

            % logicals represented as uint8
            if strcmp(char(EType), 'logical'),      EType = 'uint8';        end

            fwrite(fid,dataAll,char(EType));
            fclose(fid);

            f=fid;

        end
        
    end
    
end
