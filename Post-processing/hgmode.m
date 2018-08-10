% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% 
% 10/22/10 : Correct a missing scale factor of a0px, a0py etc.
% 05/18/10 : calculate the HG mode at z=0 point for various orders
% on a grid X,Y. 
% NB: This code will fail for m, n > 80, due to finite precision.
% 05/19/10 : Correct for exponents to be separated as a0x != a0y generally.
function mode_XY = hgmode(a0x,a0y,m,n,X,Y)
	Y = reshape(Y,size(X))';
	[XX,YY]=meshgrid(X,Y);
	scale = exp(0.5*(m*log(2)+n*log(2)+gammaln(m+1)+gammaln(n+1) + 2*0.5*log(pi) + log(a0x) + log(a0y)));
	mode_XY = exp(-(XX.^2/(2*a0x^2))-(YY.^2/(2*a0y^2))).*kron(hermitepoly(m,X/a0x)/sqrt(scale),hermitepoly(n,Y/a0y)/sqrt(scale));
	return
