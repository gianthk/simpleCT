function Y=imresizen(X,scaling,varargin)
%Resize an n-dimensional array
%
%      Y=imresizen(X,scaling,extrapMethod)
%
%IN:          
%             
%          X: n-dimensional Input array
%    scaling: scaling factor(s) as scalar or n-vector. Scalings along
%             singleton dimensions of X are always ignored.
%   extrapMethod: extrapolation method (same options as for
%                 griddedInterpolant). Default is linear interpolation.
%
%OUT:         
%             
%    Y: resized array
%
%EXAMPLES:
%
%       img4D=rand(30,20,10,40);
%
%    (1)   imresizen(img4D,2);
%    (2)   imresizen(img4D,1./[3,2,1,4]);
%    (3)   imresizen(___, 'cubic')


if ~isa(X,'single') && ~isa(X,'double')
    error 'Input array must be single or double'
end

N=ndims(X);

scaling(1,1:N)=scaling(:).';

sz=size(X);

xvec=cell(1,N);
yvec=cell(1,N);
szy=nan(1,N);
nonsing=true(1,N);

for i=1:N
 
     n=sz(i);
     
     if n==1 %for vector input
         
         nonsing(i)=0;
         szy(i)=1;
         continue
     end
     
     szy(i)=round(sz(i)*scaling(i));  
      m=szy(i);
     
     xax=linspace(1/n/2, 1-1/n/2 ,n); 
       xax=xax-.5;
       
     yax=linspace(1/m/2, 1-1/m/2 ,m); 
       yax=yax-.5;
   
   xvec{i}=xax; 
   yvec{i}=yax;     
    
end

 xvec=xvec(nonsing);
 yvec=yvec(nonsing);

 F=griddedInterpolant(xvec,squeeze(X),varargin{:});
 Y=reshape(F(yvec),szy);