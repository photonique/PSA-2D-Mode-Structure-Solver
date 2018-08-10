% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%  
%  11/23/11: Use the refractive index in k=(2 Pi n)/lambda
%  
%  11/17/11: Include the overlap integral, with a HGm,n beam offset from
%            the z=0 beam waist. We believe this is could be imporant.
% 
%  04/19/10: Calculate the overlap integral with scale gaussian beam.
%  05/18/10: Allow to find overlap with higher order modes. We generalize
%            the overlap calculation for arbitrary higher order modes m,n.
%            Caveat : for very large mode number m,n > 80 numerical
%            overflow occurs.
% 05/19/10:  The signal x,y spot sizes are not required and commented out
%            as they are potentially miselading.
% 
% 11/16/10:  Force the eigenmode & HG-mode to be unity energy / orthonormal
% 
% 12/01/10:  Right way of normalizing grid points is by using step size =
% range / (Npts - 1), as Npts on grid form only Npts-1, sections.
% Computational science 101!
% 
function overlap = spotsize_scalefactor_woffset(a0px,a0py,mode, m, n, zoffset, lambda,nref)
if ( nargin < 3 )
    error('Usage: overlap = spotsize_scalefactor(a0px,a0py,mode,m=0,n=0,zoffset=0,lambda=1560e-9,n=1.78)');
end
if ( nargin < 4 )
    m = 0; n = 0;
end
if ( nargin < 6 )
    zoffset = 0; 
    lambda=1560e-9;
    nref=1.0;
end

% for older data @ high-resolution prior to Aug 2010
%Xpts = 513; 1024; Yplambda=1560e-9;ts = Xpts; x = linspace(-0.8e-3, +0.8e-3,Xpts); y = x;
%dx = 1.6e-3/(Xpts-1); dy = 1.6e-3/(Ypts-1);

% for XXXX data sets of Oct-Nov 2011 as well.
% for older modeling upto September 2010
Xpts = 501; Ypts = Xpts; x = linspace(-1.6e-3,+1.6e-3,Xpts); y = x;
dx = 3.2e-3/Xpts; dy = 3.2e-3/Ypts;

% for modeling NU-experiment data since October 2010
%Xpts = 512; Ypts = Xpts; x = linspace(-371.0938e-6/2, +371.0938e-6/2,Xpts); y = x;
%dx = 371.0938e-6/Xpts; dy = 371.0938e-6/Ypts;

%a0px = 200e-6/sqrt(10); a0py = 100e-6/sqrt(10);
%a0sx = sqrt(2)*a0px; a0sy = sqrt(2)*a0py;

% exact fundamental order HG-modes!
% if ( m + n < 1 )
%  [XX,YY] = meshgrid(x,y);
%  f1 = 1/(sqrt(pi*a0px*a0py));
%  G0 = f1*exp(-XX.^2/(2*a0px^2)).*exp(-YY.^2/(2*a0py^2));
% else
%if ( m == 0 && n == 0 )  
%  G0 = exp(-x.^2/(2*a0px^2))'*exp(-y.^2/(2*a0py^2))./sqrt(pi*a0px*a0py);  
%
G0 = genhgmode(a0px,a0py,m,n,x,y,zoffset,lambda,nref);

% end
%  mesh(x,y,abs(G0).^2);
%  view(0,90)
%  figure
%  mesh(x,y,abs(mode).^2);
%  view(0,90)
s1 = sum(sum(abs(G0).^2))*dx*dy;
s2 = sum(sum(abs(mode).^2))*dx*dy;
%[s1, s2]
%assert( s1 > 0.75 && s2 > 0.9 )
%overlap = sum(sum(abs(mode.*G0).^2))*dx*dy/(s1*s2);
overlap = abs(sum(sum(mode.*conj(G0)))*dx*dy).^2/(s1*s2);
return
