% 12/07/09 :  (C) 2009 Muthiah Annamalai
%             Nonlinear Optics & Nanophotonics Lab,
%             University of Texas at Arlington.
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% get odd-and even block matrices of size My x My 
% embedded in the chequered matrix C of size Mx x Mx.
% No symmetry assumptions made on C.
% see also: blockset, chqmat
function [Ce,Co]=blockget( C, Mx, My )
   Ce = zeros(My*ceil(Mx/2)); Co = zeros(My*floor(Mx/2));
   rows_C = size(C,1);
   ri = 1;  ci = 1;
   for idx = 1:2:(rows_C/My)+ mod(rows_C/My,2)
       ci = 1;
       for idy = 1:2:(rows_C/My)+ mod(rows_C/My,2)

           IS1 = 1 + (idx-1)*My:idx*My; IS2 = 1 + (idy-1)*My:idy*My;
           IE1 = 1 + (ri-1)*My:ri*My; IE2 = 1 + (ci-1)*My:ci*My;
           Ce( IE1, IE2 ) = C( IS1, IS2 );
           
           if ( max(idx + 1,idy+1) <= (rows_C/My) )
             IS1 = 1 + (idx+1-1)*My:(1+idx)*My; IS2 = 1 + (idy+1-1)*My:(1+idy)*My;
             IE1 = 1 + (ri-1)*My:ri*My; IE2 = 1 + (ci-1)*My:ci*My;
             Co( IE1, IE2 ) = C( IS1, IS2 );
           end
           
           ci = ci + 1;
       end
       ri = ri + 1;
   end

%    ri = 1;  ci = 1;
%    for idx = 2:2:(rows_C/My)
%        ci = 1;
%        for idy = 2:2:(rows_C/My)
% 
%            IS1 = 1 + (idx-1)*My:idx*My; IS2 = 1 + (idy-1)*My:idy*My;
%            IE1 = 1 + (ri-1)*My:ri*My; IE2 = 1 + (ci-1)*My:ci*My;
%            Co( IE1, IE2 ) = C( IS1, IS2 );
%            
%            ci = ci + 1;
%        end
%        ri = ri + 1;
%    end
   
   return
   %end

%   rows_C = size(C,1); cols_C = size(C,2);
%  Ce = zeros(ceil(rows_C/My/2)*My); Co = zeros(floor(rows_C/My/2)*My);
%    for i = 1:2: ( ceil(rows_C/My/2) + mod(rows_C/My,2) )
%     for j = 1:2: ( ceil(rows_C/My/2) + mod(rows_C/My,2) )
%       % extract even position
%       IS1 = 1+(ceil(j/2)-1)*My:ceil(j/2)*My; IS2 = 1+(ceil(i/2)-1)*My:ceil(i/2)*My;
%       IE1 = 1+(i-1)*My:i*My; IE2 = 1 +(j-1)*My:j*My;
%       Ce(IS1,IS2)   = C(IE1,IE2);
%       Ce(1+(ceil(j/2)-1)*My:ceil(j/2)*My, 1+(ceil(i/2)-1)*My:ceil(i/2)*My)  = C(1 +(j-1)*My:j*My,1+(i-1)*My:i*My);
%       
%       % odd position
%       i2 = i + 1; j2 = j + 1;
%       if ( max(i2 ,j2) <= floor(rows_C/My/2) + mod(rows_C/My,2) )
%         Co(1+(i2/2-1)*My:i2/2*My,1+(j2/2-1)*My:(j2/2)*My) = C(1 +(i2-1)*My:i2*My,1+(j2-1)*My:j2*My);
%         Co(1+(j2/2-1)*My:(j2/2)*My,1+(i2/2-1)*My:i2/2*My) = C(1 +(i2-1)*My:i2*My,1+(j2-1)*My:j2*My);
%       end
%       
%     end
%    end
