%% 
%% (C) 2010 Muthiah Annamalai <muthuspost@gmail.com>
%%     Nonlinear Optics & Nanophotonics Lab, UT-Arlington.
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%% 
%% 
%% 08/03/10 : Code for visualize the eigenmodes of higher order pump 
%%            with pump power kW for several polynomials over a
%%            grid size of 3.2mm x 3.2mm square area for 2cm PSA, and at
%%            z=0 distance (similar to Apr/May '10 Eigenmode data
%%            calculations of the MV9 solver done at UTA).
%%            
%%            The grid information is sampled at 501 points centered
%%            around 0-point; ofcourse this data was generated from
%%            Mathematica for signal mode with 
%%            waist size (sqrt(2) over pump) & higher order HG-basis.
%% 
%% 08/11/10 : Reorganize code to delineate where Mx,My,a0p{x,y}
%%            correspond to which data files. Eigenmodes are generated well!
%% 
%% 08/27/10 : Looks like from test_ortho.m, no red-flags in the spatial 
%%            eigenmodes for several higher order pump cases, which are
%%            very orthogonal.
%% 
clear all
close all
more off

addpath('BBN-HG-basis-modes-data/')

if ( 0 )
  % problem parameters : 2cm PSA crystal
  Mx = 128; My = 128; a0px = 200e-6; a0py = 200e-6; Lz =  20e-3; P0 = 16.25e3;

  % HG-basis functions are required on the z=0cm 
  % PSA o/p plane.
  hermiX = load('hermiX513_200.dat');
  hermiY = load('hermiX513_200.dat');
  
  q= 4; m = 4;
  
  % 200x200 spot size @ 256x64 modes for 17.25kW power
  % load 'MV11_Eigenvalue_q00_m00_Mx128_My128-P016250-A0x200-A0y200-22-Jul-2010.mat'
  % load 'MV11_Eigenvalue_q02_m02_Mx128_My128-P016250-A0x200-A0y200-24-Jul-2010.mat'
  load 'MV11_Eigenvalue_q04_m04_Mx128_My128-P016250-A0x200-A0y200-28-Jul-2010.mat'
end

if ( 0 )
  % for the 400x100 spot size
  Mx = 256; My = 64; a0px = 400e-6; a0py = 100e-6; Lz = 20e-3; P0 = 17.25e3;

  hermiX = load('hermiX513_400.dat');
  hermiY = load('hermiX513_100.dat');

  q= 6; m = 0;

  % 400x100 spot size @ 256x64 modes for 17.25kW power
  %load 'MV11_Eigenvalue_q00_m00_Mx256_My64-P017250-A0x400-A0y100-04-Aug-2010.mat';
  %load 'MV11_Eigenvalue_q04_m00_Mx256_My64-P017250-A0x400-A0y100-24-Jul-2010.mat'
  load 'MV11_Eigenvalue_q06_m00_Mx256_My64-P017250-A0x400-A0y100-28-Jul-2010.mat'
end


  % for the 400x100 spot size
  Mx = 32; My = 32; a0px = 100e-6; a0py = 25e-6; Lz = 20e-3; P0 = 2.4e3;

  hermiX = load('hermiX513_100.dat');
  hermiY = load('hermiX513_25.dat');

  q= 2; m = 0;

  % 200x100 spot size @ 32x32 modes for 2.4kW power
  load 'MV11_Eigenvalue_q02_m00_Mx32_My32-P02400-A0x100-A0y25-11-Aug-2010.mat'

% eigenmode are calculated in grid size of area 1.6mm x 1.6mm
Xpts = 501; Ypts = Xpts; x = linspace(-1.6e-3,+1.6e-3,Xpts); y = x;
nmodes = 25; 5;

a0sx = sqrt(2)*a0px; a0sy = sqrt(2)*a0py;

psastr = sprintf('pump-q%02d-m%02d-P0=%g-kW-spot-%dx%d-Mx%d-My%d',q,m,P0/1e3,round(a0px/1e-6), round(a0py/1e-6),Mx,My);

fprintf('Working on PSA problem:\n %s problem\n',psastr)
l=flipud(l);
evec = fliplr(evec);
vecs = evec(1:Mx*My,1:nmodes) + 1i*evec(1+Mx*My:end,1:nmodes);

eigenmodes = zeros( Ypts, Xpts, nmodes );

%% visualize the eigenvector/mode
for idmode = 1:nmodes;
     fprintf('Working on Eigenmode %g\n',idmode);
     mode = zeros(size(hermiY,2),size(hermiX,2));
     idx = 0;

    mode = [hermiX(1:Mx,:).'*reshape(vecs(:,idmode),My,Mx).'*hermiY(1:My,:)].';


   %  %% slow way to generate data, but we are sure here I believe.
if ( 0 )
     for mx = 1:Mx
         for my = 1:My
             idx = idx + 1;
             mode = mode + vecs(idx,idmode)*(hermiY(my,:).')*(hermiX(mx,:));
         end
     end
end
     %% all modes are orthonormal over the whole space, but in the
     %% limited grid size of 1.6x1.6mm^2 we make them orthogonal
     %% again.
     mode = mode/sqrt((sum(sum(abs(mode.^2)))*(1.6e-3/513).^2));
     sum(sum(abs(mode.^2)))*(1.6e-3/513).^2

     eigenmodes(:,:,idmode) = mode;

    if ( 1 )
     %figure
     surf(x(1:2:end)/1e-3,y(1:2:end)/1e-3,abs(mode(1:2:end,1:2:end)).^2)
     view(0,90)
     axis square
     colorbar
     shading interp
     title(sprintf('Mode# %d, [Eigval=%g] , %s\\mum^2',idmode-1,l(idmode),psastr))
     xlabel(['x axis (mm),', num2str(Mx),'HG modes']); 
     ylabel(['y axis (mm),', num2str(My),'HG modes']);
     set(gcf,'Color','White')
     print(gcf,'-dpng',['HO-Spatial-Eigenmode-',num2str(idmode-1),'-Eigenvalue=',num2str(l(idmode)),psastr,'.png']);
    end
end
save('-mat',['MV11_Eigenmodes_',psastr,'.mat'],'eigenmodes','l','evec','x')
fprintf('%s\n',['Saved file : MV11_Eigenmodes_',psastr,'.mat'])
