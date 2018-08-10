% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% 11/23/11 : Use refractive index argument
% 11/17/11 : Like kvhgmode2d - but with elliptic spot size allowed.
% 10/22/10 : Correct a missing scale factor of a0px, a0py etc.
% 05/18/10 : calculate the HG mode at z=0 point for various orders
% on a grid X,Y. 
% NB: This code will fail for m, n > 80, due to finite precision.
% 05/19/10 : Correct for exponents to be separated as a0x != a0y generally.
function mode_XY = genhgmode(a0x,a0y,m,n,X,Y, z, lambda, nref)
    if ( nargin < 9 )
	nref = 1.0;
    end

    zrx = 2*pi*nref*a0x^2/lambda;
    zry = 2*pi*nref*a0y^2/lambda;
    
    azx = a0x*sqrt(1 + (z/zrx)^2 );
    azy = a0y*sqrt(1 + (z/zry)^2 );    
    phizx = atan( z/zrx ); phizy = atan( z/zry );    
    k = 2*pi*nref/lambda;
    
    Y = reshape(Y,size(X))';
	[XX,YY]=meshgrid(X,Y);
    
    if ( abs( z ) > 5*eps )
        rzx = z*(1 + (zrx/z)^2);
        rzy = z*(1 + (zry/z)^2);
        scaleExp = exp(-i*(m+0.5)*phizx -i*(n+0.5)*phizy);        
        scaleExp = scaleExp*exp(+i*k*(((XX.^2)./(2*rzx))+((YY.^2)/(2*rzy))));
    else
        scaleExp = 1;
    end   
    
	scale = exp(0.5*(m*log(2)+n*log(2)+gammaln(m+1)+gammaln(n+1) + 2*0.5*log(pi) + log(azx) + log(azy)));
	mode_XY = exp(-(XX.^2/(2*azx^2))-(YY.^2/(2*azy^2))).*kron(hermitepoly(m,X/azx)/sqrt(scale),hermitepoly(n,Y/azy)/sqrt(scale));    
    mode_XY = mode_XY.*scaleExp;
	return
