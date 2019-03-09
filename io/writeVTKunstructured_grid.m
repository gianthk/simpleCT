function [ check ] = writeVTKunstructured_grid( filename, points, point_data )
%WRITEVTKUNSTRUCTURED_GRID writes unstructured_grid VTK file..
%
%   INPUTS:     points (n x 3)      data points coordinates
%       
%   ______________________________________________________
%
%   Author:         Gianluca Iori (gianluca.iori@charite.de)
%   BSRT - Charite Berlin
%   Created on:   14/03/2017
%   Last update:  14/03/2017
%   ______________________________________________________
%

    check = false;
    
    if nargin < 3,      point_data = [];       end
    
    if size(points,2) ~= 3
        error('points array must be n x 3..');
    end
    
    if ~isempty(point_data)
        if size(point_data,2) ~= 3
            error('point_data array must be n x 3..');
        end
        if size(point_data,1) ~= size(points,1)
            error('point_data array must have the same number of rows as points..');
        end
        
    end
    
    [PATHSTR,NAME,EXT] = fileparts(filename);
    
    % open VTK file
    fid = fopen(filename, 'w+');
    
    % header
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Unstructured Grid %s\n', NAME);
    fprintf(fid, 'ASCII\n');
    fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
    
    % points
    fprintf(fid, 'POINTS %i float\n', size(points,1));
    for i=1:size(points,1)
        fprintf(fid, '%f %f %f\n', points(i,1), points(i,2), points(i,3));
    end
    fprintf(fid, '\n');
    
    % cells
    fprintf(fid, 'CELLS %i %i\n', (size(points,1)-1), (size(points,1)-1)*3);
    for i=1:(size(points,1)-1)
        fprintf(fid, '2 %i %i\n', i-1, i);
    end
    fprintf(fid, '\n');
    
    % cell_types
    fprintf(fid, 'CELL_TYPES %i\n', (size(points,1)-1));
    for i=1:(size(points,1)-1)
        fprintf(fid, '3\n');
    end
    fprintf(fid, '\n');
    
    if ~isempty(point_data)
        % point_data
        fprintf(fid, 'POINT_DATA %i\n', size(points,1));
        fprintf(fid, 'VECTORS displacement float\n');
        for i=1:size(points,1)
            fprintf(fid, '%f %f %f\n', point_data(i,1), point_data(i,2), point_data(i,3));
        end
        fprintf(fid, '\n');
    end
    
    fclose(fid);
    
end

