%% 
%% (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
%%     Nonlinear Optics & Nanophotonics Lab, UT-Arlington.
%%
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
%
%% Usage: c = kronsum( a, b )
%%
%% Implements kronecker sum of 'a', 'a' ala direct sum
%% just like the kronecker product operator ordering.
%%
%% The second matrix/vector 'b', is added to every element in the 
%% first matrix/vector 'a', resulting in the combined matrix 'c'.
%% 
%% eg: kronsum( [ 0 1 2 ]' , [ 0 1 2] ) =
%% 
%%              [ 0 1 2; 
%%                1 2 3;
%%                2 3 4 ]                
%%
function c = kronsum( a, b )
  c = kron( ones(size(b)), a ) + kron( b, ones(size(a)) );
  return
end
