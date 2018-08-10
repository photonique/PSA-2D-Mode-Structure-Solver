% 12/07/09 :  (C) 2009 Muthiah Annamalai
%             Nonlinear Optics & Nanophotonics Lab,
%             University of Texas at Arlington.
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
%
% set odd-and even block matrices of size My x My 
% embedded into the chequered matrix C of size Mx x Mx.
% No symmetry assumptions made on C.
% see also: blockget, chqmat
function C=blockset( Ce,Co, Mx, My )
   C = zeros(Mx*My);   
   rows_C = size(C,1);
   ri = 1;  ci = 1;
   for idx = 1:2:(rows_C/My)+ mod(rows_C/My,2)
       ci = 1;
       for idy = 1:2:(rows_C/My)+ mod(rows_C/My,2)

           IS1 = 1 + (idx-1)*My:idx*My; IS2 = 1 + (idy-1)*My:idy*My;
           IE1 = 1 + (ri-1)*My:ri*My; IE2 = 1 + (ci-1)*My:ci*My;
            C( IS1, IS2 ) = Ce( IE1, IE2 );
           
           if ( max(idx + 1,idy+1) <= (rows_C/My) )
             IS1 = 1 + (idx+1-1)*My:(1+idx)*My; IS2 = 1 + (idy+1-1)*My:(1+idy)*My;
             IE1 = 1 + (ri-1)*My:ri*My; IE2 = 1 + (ci-1)*My:ci*My;
             C( IS1, IS2 ) = Co( IE1, IE2 );
           end           
           ci = ci + 1;
       end
       ri = ri + 1;
   end
   return
   %end
   
