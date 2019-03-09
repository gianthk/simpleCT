classdef medtooltable < handle
    %MEDTOOLTABLE class for imput and output of a table formatted for the use with
    %   Medtool from Dr.Pahr Ingenieurs e.U.
    %
    %   see:    http://www.dr-pahr.at/software_en.php
    %
    %   static methods are implemented to read and write the table in the following formats:
    %   -   Comma Separated Values (.CSV)
    %   -   Text (.txt)
    %
    %   The property TABLE is a MATLAB structure with the fields:
    %       header      (1 x n) Cell Array containing the names of the
    %                   table columns. These are stored in string format
    %                   and should start with "$"
    %       
    %       "field"     one table field for each element of the cell array "header"
    %                   with the name specified by the header itself
    %
    %   EXAMPLE:
    %           tab = medtooltable.csvreadtable
    %       or:
    %           tab = medtooltable;
    %           tab.readtable;
    %           tab.load;
    %   ______________________________________________________
    %
    %   Author:         Gianluca Iori (gianthk.iori@gmail.com)
    %   BSRT - Charite Berlin
    %   Created on:     04/08/2017
    %   Last update:    10/12/2017
    %
    %   this class is part of the synchro toolbox
    %   ______________________________________________________
    
    properties (Access = public)

        table;              % MATLAB structure with field: "header" and table columns with names as specified in "header"
        
        % delimiter;          % ',' ';'(default) ' '{space} '.' ':'
        
    end
    
    methods (Static = true)
        
        function [table, pathname] = readtable(filename)
            %READTABLE reads file containing a medtool table.
            %   The first row is assumed to contain table headers.
            %   tab = readtable                     opens file selection gui and reads table
            %   tab = readtable(filename)           reads table data
            %   ______________________________________________________
            %
            %   Author:         Gianluca Iori (gianluca.iori@charite.de)
            %   BSRT - Charite Berlin
            %   Created on:   27/06/2016
            %   Last update:  07/12/2017
            %   ______________________________________________________
            %
            %   See also TXTREAD
            %   ______________________________________________________
            
            if ~exist('filename','var') || isempty(filename)
                [filename, pathname] = uigetfile({ '*.csv;*.CSV;*.txt;*.TXT','Supported Files (*.csv,*.CSV,*.txt,*.TXT)';...
                                                   '*.csv;*.CSV','Comma Separated Values (*.csv,*.CSV)';...
                                                   '*.txt;*.TXT','Text files (*.txt,*.TXT)';...
                                                   '*.*', 'All Files (*.*)'},'Load Medtool table','MultiSelect','off');
                                               
                if filename==0
                    % user pressed cancel
                    return;
                end
               
                filename = [pathname filename];
            else
                % get input fileparts
                [pathname,name,ext] = fileparts(filename);
            end

            fid = fopen(filename);

            % get table headers and table width
            firstline_ = fgets(fid);
            firstline = regexp(firstline_,'[;,]','split');
            firstline{length(firstline)} = firstline{length(firstline)}(1:end-2);
            table.header = firstline;
            column = length(firstline);

            % get columns format
            line_ = fgets(fid);
            line = regexp(line_,'[;,]','split');
            format = '';                        % initialize format string for textscan
            for i = 1:length(line)
%                 if ~all(isstrprop(strrep(line{i}, '.', ''), 'digit'))
%                 if isempty(str2num(strrep(line{i}, '-', '_')))
                if isempty(str2num(line{i}))
                    format = [format '%s'];     % string read
                elseif round(str2num(line{i})) == str2num(line{i})
                    format = [format '%d'];     % syntax for reading integer
                else
                    format = [format '%f'];     % float read
                end
            end

            % read table contents
            fseek(fid,0,-1);                    % return to beginning of file
            fgetl(fid);                         % skip header line
            tabledata = textscan(fid,format,'delimiter',{',',';'});       % textscan table

            fclose(fid);

            % assemble output struct
            for i=1:length(line)
                if isempty(firstline{i}),   continue;   end
                table.(regexprep(strrep(firstline{i},'$',''),'\s','')) = tabledata{i};
            end
        end
        
        function writetable(table, filename)
            %WRITETABLE writes table to comma separated file.
            %   The first row is assumed to contain table headers.
            %   writetable(table)                opens output filename selection gui and writes table
            %   writetable(table, filename)      writes table to given file
            %   ______________________________________________________
            %
            %   Author:         Gianluca Iori (gianluca.iori@charite.de)
            %   BSRT - Charite Berlin
            %   Created on:   27/12/2016
            %   Last update:  04/08/2017
            %   ______________________________________________________
            
            %% check inputs
            if nargin < 1
                error('Not enough inputs.   Example call: txtwritetable(table, filename)');
            end

            if nargin < 2,      filename='';      end

            %% select output file
            if ~exist('filename','var') || isempty(filename)
                % output file
                [filename, pathname] = uiputfile({ '*.csv;*.CSV;*.txt;*.TXT','Supported Files (*.csv,*.CSV,*.txt,*.TXT)';...
                                                   '*.csv;*.CSV','Comma Separated Values (*.csv,*.CSV)';...
                                                   '*.txt;*.TXT','Text files (*.txt,*.TXT)';...
                                                   '*.*', 'All Files (*.*)'},'Load Medtool table','MultiSelect','off');
                                               
                if filename==0
                    % user pressed cancel
                    return;
                end
               
                filename = [pathname filename];
                
            end

            %% open output TXT file
            fid = fopen(filename, 'w');

            %% write header line
            % print header line
            headerstr = table.header{1};
            for i = 2:length(table.header)
                headerstr = strcat(headerstr, ';', table.header{i});
            end
            % headerstr = strcat(headerstr, '\n');

            % fprintf to txt file
            fprintf(fid,'%s\n',headerstr);

            %% write table content
            for row = 1:length(table.(strrep(table.header{1},'$','')))
                % print line content
                if iscell(table.(strrep(table.header{1},'$','')))
                    linestr = table.(strrep(table.header{1},'$','')){row};
                else
                    linestr = num2str(table.(strrep(table.header{1},'$',''))(row));
                end

                for i = 2:length(table.header)
                    if isempty(table.header{i}),   continue;   end
                    if iscell(table.(strrep(table.header{i},'$','')))
                        linestr = strcat(linestr, ';', table.(strrep(table.header{i},'$','')){row});
                    else
                        linestr = strcat(linestr, ';', num2str(table.(strrep(table.header{i},'$',''))(row)));
                    end
                end
                fprintf(fid,'%s\n',linestr);

                clear linestr;

            end

            %% close file
            fclose(fid);
            fprintf('table data written to file.\n');
            
        end
        
        function [header] = getHeader(table)
            if isfield(table, 'header')                 % header is already in table
                warning('medtooltable.getHeader: table contains header field already!');
                return;
            else
                header = fieldnames(table)';        % get field names
                for i = 1:length(header)
                    header(i) = strcat('$', header(i));
                end
            end
        end
        
    end
    
    methods
        % input/output
        function load(this, filename)
            % function load(this, filename)
            % load medtooltable
            %
            % supported file formats:    - Comma Separated Values   (*.csv; *.CSV)
            %                            - Text files               (*.txt,*.TXT)
            %
            % example call:     tab = sct1.io.medtooltable
            %                   tab.load

           % get input file if no argument is given
           if nargin < 2
               
               [filename, pathname] = uigetfile({  '*.csv;*.CSV;*.txt;*.TXT','Supported Files (*.csv,*.CSV,*.txt,*.TXT)';...
                                                   '*.csv;*.CSV','Comma Separated Values (*.csv,*.CSV)';...
                                                   '*.txt;*.TXT','Text files (*.txt,*.TXT)';...
                                                   '*.*', 'All Files (*.*)'},'Load Medtool table','MultiSelect','off');
               if filename==0
                    % user pressed cancel
                    return;
               end
               
               if isempty(filename)
                   if length(filename) > 1
                       if ~ischar(filename{1})
                           this = [];
                           fprintf('opening aborted');
                           msg  = 'opening aborted';
                           return;
                       end
                   elseif ~ischar(filename)
                       this = [];
                       fprintf('opening aborted');
                       msg  = 'opening aborted';
                       return;
                   end
               end
               
               filename = [pathname filename];
           end
           
           this.table = this.readtable(filename);
           
        end
    
        % getters
        
        % setters
        function setHeader(this)
            this.table.header = medtooltable.getHeader(this.table);
        end
        
        function setTable(this, table)
            this.table = table;
            this.setHeader;
        end
        
        % medtooltable class constructor
        function this = medtooltable()
            this.table = []; 
        end
        
    end
        
end