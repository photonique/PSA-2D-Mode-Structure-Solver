% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% 07/13/10 : Calculate the Lambda_{m,l,n}(z=0) for higher order pump
%            Elementary method works very well, but limited to My <= 80.
%            Superior method using the gammaln is buggy which is bad,
%            because if it worked we can explore arbitrary size modes My.

% 
% m : pump mode
% My: pump mode along Y
% a0y: pump spot size along Y
% 

function overlap = calc_higherorder_overlap( m, My )
for l = 0:My-1
 for n = 0:My-1
  scale = exp( 0.5*(gammaln(l+1) + gammaln(n+1) ...
               - gammaln( m  + 1) - (m+1)*log(2) - 7/2*log(pi)));

  val  = 0;
  % selection rule + traditional way of calculating the integrals.
  if ( mod( m + l + n , 2 ) == 0 )
   for v = 0:min(l,n)
    w = (m + n + l - 2*v + 1)/2;
    for k = 0:m
     val = val + nchoosek( m, k)*gamma(w-(m-k))/gamma(v+1)*gamma(w-k)/gamma(l-v+1)*gamma(w-(l+n-2*v))/gamma(n-v+1);
    end
   end
  end
  val = val*scale;

  overlap(l+1,n+1) = val;
 end
end

return
end
