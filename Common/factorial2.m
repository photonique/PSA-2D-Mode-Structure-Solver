%% 
%% (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
%%     Nonlinear Optics & Nanophotonics Lab, UT-Arlington.
%% 
%% This code may be used or distributed under terms of MIT License.
%% This file is part of the PSA-2D-Mode-Structure-Solver project.
%%
%% Fast evaluation of double-factorial aka 'Factorial2' function following
%% definitions from Arfken, Wolfram-Mathworld for odd-case.
%% 
%% Calculate the factorial2 function for a matrix/vector 
%% input n. Currently only dimension condense / calculate /expand
%% procedure works fast.
%% 
%% Use gammaln functions but this is useless for large values 
%% as is the factorial() commands. So just report the error.
%% 
function v = factorial2( n )
    if ( max(max(n)) > 300 )
        error('Factorial2 function overflows for n > 300')
    end
    
    s = size( n );	
    n = reshape(n,1,numel(n));
    v = 0*n;
    
    [tmp,id0]=find( n==0);
    v(id0) = 1;
    clear id0, tmp;
    
    [tmp,idm1] = find( n==-1);
    v(idm1) = 1;
    clear tmp, idm1;    
    
    [tmp,idod] = find( mod(n,2) == 1 );
    %%v(idod) = gamma( n(idod)/2 + 1).*2.^((n(idod)+1)/2)/sqrt(pi);
    v(idod) = exp(gammaln(n(idod)/2 + 1) + (n(idod)+1)/2*log(2) - 0.5*log(pi));
    clear tmp, idod;
    
    [tmp,idev] = find( mod(n,2) == 0 );
    %v(idev) = 2.^(n(idev)/2).*gamma( 1 + n(idev)/2);
    v(idev) = exp( n(idev)/2*log(2) + gammaln(1 + n(idev)/2) );
    clear tmp, idev;
			
    v = reshape( v, s(1),s(2));
    return
end
