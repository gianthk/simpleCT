function [DICOMslicenames] = getDICOMslicelist(DICOMfile, slice_in, slice_end)
%GETDICOMSLICELIST creates list of DICOM files given start and end slices
%   list = getDICOMslicelist(DICOMfile, slice_in, slice_end) create list of
%   DICOM slices given a filename and initial and final slices
%   ______________________________________________________
%
%   Author:         Gianluca Iori <gianthk.iori@gmail.com>
%   Created on:     22/02/2018
%   Last update:    22/02/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    count = 1;
    DICOMfile = regexp(DICOMfile, '\.', 'split');	% split filename and file format string
    formatstring = ['.' DICOMfile{2}];              % file format string
    DICOMfile = DICOMfile{1};                       % file name (without format string)
    for slice = [slice_in : slice_end]
        ndigits = numel(num2str(slice));
        DICOMslicenames{count} = [DICOMfile(1:end-ndigits) num2str(slice) formatstring];
        count = count + 1;
    end
end
