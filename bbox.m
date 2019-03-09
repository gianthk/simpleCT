function [row0, rowd, col0, cold, slice0, sliced] = bbox(bw, offset, dsize, graphics)
%BBOX calculates bounding box limits of input binary image 
%   [row0, rowd, col0, cold, slice0, sliced] = bbox(bw) calculates origin [row0 col0 slice0] and
%   size [rowd cold sliced] of the bounding box inscribing ones in input binary data bw
%
%   [row0, rowd, col0, cold, slice0, sliced] = bbox(bw, offset) adds offset pixels to the bbox limits
%
%   [row0, rowd, col0, cold, slice0, sliced] = bbox(bw, offset, dsize) performs imclose with disk
%   structuring element with radius dsize before calculating the bbox
%
%   [row0, rowd, col0, cold, slice0, sliced] = bbox(bw, offset, dsize, graphics) display graphical output
%
%   BBOX coordinates are returned in pixels.
%   Input bw can be 2D or 3D matrix of reals, integers or logicals but must be binary.
%   If input bw is 2D matrix BBOX will return only row and col bbox indexes.
%   ______________________________________________________
%
%   Authors:        Johannes Schneider
%                   Gianluca Iori (gianthk.iori@gmail.com)
%   BSRT - Charite Berlin
%   Created on:   26/01/2018
%   Last update:  08/02/2018
%
%   this function is part of the synchro toolbox    
%   ______________________________________________________

    if nargin < 4,                  graphics = false;           end
    if nargin < 3,                  dsize = 0;                  end
    if isempty(dsize),              dsize = 0;                  end
    if nargin < 2,                  offset = 0;                 end
    if isempty(offset),             offset = 0;                 end
    
    narginchk(1,4);
    
    validateattributes(graphics,{'numeric' 'logical'},{'nonempty','binary'},mfilename,'graphics',4);
    validateattributes(dsize,{'numeric'},{'nonempty','finite','integer'},mfilename,'dsize',3);
    validateattributes(offset,{'numeric'},{'nonempty','finite','integer'},mfilename,'offset',2);
    validateattributes(bw,{'numeric' 'logical'},{'nonempty','binary'},mfilename,'bw',1);
    
    %% remove artefacts > erode/dilate
    if dsize>0
        fprintf('bbox: use imopen with caution!\n');
        se  = strel( 'disk' , dsize );
        bw = imopen( bw , se );
        
        % imclose might cause edge structures to be reduced.. we therefore
        % need to consider a minimum offset if we perform imopen!
        if offset<=dsize
            offset = dsize+1;
        end
    end
    
    %% find bbox
    switch ndims(bw)
        case 3
            % project along each dimension
            maxXY       = squeeze(max( bw , [] , 3 ));
            maxROW      = max( maxXY , [] , 2 );
            maxCOL      = max( maxXY , [] , 1 );
            maxXZ       = squeeze(max( bw , [] , 2 ));
            maxSLICE    = max( maxXZ , [] , 1 );

            % find indexes of first and last non-zero elements
            row0 = find( maxROW > 0 , 1 , 'first' );
            row1 = find( maxROW > 0 , 1 , 'last' );

            col0 = find( maxCOL > 0 , 1 , 'first' );
            col1 = find( maxCOL > 0 , 1 , 'last' );

            slice0 = find( maxSLICE > 0 , 1 , 'first' );
            slice1 = find( maxSLICE > 0 , 1 , 'last' );

            % add offset
            row0 = row0 - offset;           rowd = row1 - row0 + offset;
            col0 = col0 - offset;           cold = col1 - col0 + offset;
            slice0 = slice0 - offset;       sliced = slice1 - slice0 + offset;

            if offset > 0
                % check if bbox exceeds image size
                if row0<1       row0=1;     end
                if col0<1       col0=1;     end
                if slice0<1     slice0=1;   end

                bw_size = size(bw);
                if row0+rowd>bw_size(1)       rowd=bw_size(1)-row0;     end
                if col0+cold>bw_size(2)       cold=bw_size(2)-col0;     end
                if slice0+sliced>bw_size(3)   sliced=bw_size(3)-slice0; end
            end
            
            if graphics
                
                figure;
                subplot 121;
                imagesc( maxXY );
                colormap gray;
                hold on;
                plot([ col0 col0 ], [1 size( bw , 1 )] , 'r' );
                plot([ col0+cold col0+cold ], [1 size( bw , 1 )] , 'r' );
                plot([ 1 size( bw , 2 ) ], [ row0 row0 ] , 'r' );
                plot([ 1 size( bw , 2 ) ], [ row0+rowd row0+rowd ] , 'r' );
                axis image;
                
                subplot 122;
                imagesc( maxXZ );
                colormap gray;
                hold on;
                plot([1 size( bw , 3 )], [ row0 row0 ], 'r' );
                plot([1 size( bw , 3 )], [ row0+rowd row0+rowd ], 'r' );
                axis image;
                
                plot([ slice0 slice0 ], [ 1 size( bw , 1 ) ], 'r' );
                plot([ slice0+sliced slice0+sliced ], [ 1 size( bw , 1 ) ] , 'r' );
                
            end
            
        case 2
            % project along each dimension
            maxROW      = max( bw , [] , 2 );
            maxCOL      = max( bw , [] , 1 );
            
            % find indexes of first and last non-zero elements
            row0 = find( maxROW > 0 , 1 , 'first' );
            row1 = find( maxROW > 0 , 1 , 'last' );

            col0 = find( maxCOL > 0 , 1 , 'first' );
            col1 = find( maxCOL > 0 , 1 , 'last' );

            % add offset
            row0 = row0 - offset;           rowd = row1 - row0 + offset;
            col0 = col0 - offset;           cold = col1 - col0 + offset;
            slice0 = [];                    sliced = [];
            
            if offset > 0
                % check if bbox exceeds image size
                if row0<1       row0=1;     end
                if col0<1       col0=1;     end
            
                bw_size = size(bw);
                if row0+rowd>bw_size(1)       rowd=bw_size(1)-row0;     end
                if col0+cold>bw_size(2)       cold=bw_size(2)-col0;     end
            end
            
            if graphics
                
                figure;
                imagesc( bw );
                colormap gray;
                hold on;
                plot([ col0 col0 ], [1 size( bw , 1 )] , 'r' );
                plot([ col0+cold col0+cold ], [1 size( bw , 1 )] , 'r' );
                plot([ 1 size( bw , 2 ) ], [ row0 row0 ] , 'r' );
                plot([ 1 size( bw , 2 ) ], [ row0+rowd row0+rowd ] , 'r' );
                axis image;     % axis xy; 
                
            end
    end

end

