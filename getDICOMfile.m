function [DICOMfilename] = getDICOMfile(DICOMpath)
%GETDICOMFILE gets DICOM file name (without file format string)
%   [filename] = getDICOMfile(path) searches the given path for DICOM files
%   and returns one filename from the list of files contained
%   ______________________________________________________
%
%   Author:         Gianluca Iori <gianthk.iori@gmail.com>
%   Created on:     22/02/2018
%   Last update:    22/02/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    D = dir(DICOMpath);                                         % get folder contents
    DICOMfilename = D(10).name;                                 % take one file (skip first strings)
    %   DICOMfilename = regexp(DICOMfilename, '\.', 'split');	% split filename and file format string
    %   formatstring = ['.' DICOMfilename{2}];               	% file format string
    %   DICOMfilename = DICOMfilename{1};                   	% filename
end

