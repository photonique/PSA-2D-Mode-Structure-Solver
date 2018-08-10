% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%%
%% 05/17/11 : Redo same modevisualizations for data set 800x25um^2 pump spot,
%%            at several pump powers, 5kW, 10kW respectively over same grid
%%            sizes as before [08/03/10].
%%
%% 08/03/10 : Code for BBN requirements to visualize the eigenmodes of pump 400x25 spot size 
%%            with pump power {10,5,3}kW for 512x32 polynomials over a
%%            grid size of 1.6mm x 1.6mm square area for 2cm PSA, and at z=+L/2 distance.
%%            The grid information is sampled at 513 points centered
%%            around 0-point; ofcourse this data was generated from
%%            Mathematica for signal mode with waist size (sqrt(2) over pump) & higher order
%%            HG-basis.
%% 
%% 08/10/10 : See how modes reconstruct in the traditional way with the order
%%            of the loops being changed.
%%            Well looks like the fast-bilinear sums do the exact same thing
%%            as loops so we were right all along!
%% 
%% 05/25/11 : Modifications to view the eigenmodes for 800x25umsq spot
%%            and all along we are only interested in the mode-shape at
%%            the point z=0.
%% 07/18/11 : Postproc the 19.5x19.5 um^2 pump, at 185W, d_eff = 8.0pm/V
%%            for Harris Experimental parameters
%% 

clear all
close all
more off

% problem parameters : 2cm PSA crystal
% Mx = 512; 128; My = 32; a0px = 400e-6; 200e-6; a0py = 50e-6; 25e-6; Lz =  20e-3;
% Mx = 2048; My = 8; a0px = 800e-6; a0py = 25e-6; Lz =  20e-3;
Mx = 32; My = 32; a0px = 19.5e-6; a0py = a0px; Lz = 25e-3;

% eigenmode are calculated in grid size of area 1.6mm x 1.6mm
Xpts = 513; Ypts = Xpts; %x = linspace(-0.8e-3,+0.8e-3,Xpts); y = x;
%x = linspace(-1.6e-3,+1.6e-3,Xpts); y = x;
x = linspace(-0.6e-3,+0.6e-3,Xpts); y = x;
dx = 1.2e-3/(Xpts-1); dy = dx;
nmodes = 50;

a0sx = sqrt(2)*a0px; a0sy= sqrt(2)*a0py;
%addpath('MV9-EIGENVAL-DATA/')
addpath('BBN-HG-basis-modes-data/')

% HG-basis functions are required on the z=1.0cm 
% PSA o/p plane.

% if ( 0 )
% 
%   P0 = 5e3; 10e3;
%   
%   load(sprintf('BBN_z0.cm_ShermiX513_Npts513_%g.mat',round(a0px/1e-6)));
%   hermiX = Expression;
%   
%   load(sprintf('BBN_z0.cm_ShermiX513_Npts513_%g.mat',round(a0py/1e-6)));
%   hermiY = Expression;
% 
%   %load 'MV9_Eigenvalue_Mx512_My32-P010000-A0x400-A0y25-02-Aug-2010.mat'
%   %load 'MV9_Eigenvalue_Mx512_My32-P05000-A0x400-A0y25-03-Aug-2010.mat'    
%   
% %   P0 = 3e3;
% % 
% %   load('BBN_z1.cm_ShermiX513_Npts513_565.685.mat');
% %   hermiX = Expression1;
% % 
% %   load('BBN_z1.cm_ShermiX513_Npts513_35.355.mat');
% %   hermiY = Expression1;
% % load('UTA_grid1.2Mx513_Npts513_19.5.mat')
% %   %load 'MV9_Eigenvalue_Mx512_My32-P010000-A0x400-A0y25-02-Aug-2010.mat'
% %   %load 'MV9_Eigenvalue_Mx512_My32-P05000-A0x400-A0y25-03-Aug-2010.mat'
% %   load 'MV9_Eigenvalue_Mx512_My32-P03000-A0x400-A0y25-02-Aug-2010.mat'
% 
% elseif( 0 )
%   
%   P0 = 10e3;
%   
%   load('BBN_z1.cm_ShermiX513_Npts513_565.685.mat');
%   hermiX = Expression1;
% 
%   load('BBN_z1.cm_ShermiX513_Npts513_70.71.mat');
%   hermiY = Expression1;
% 
%   load 'MV9_Eigenvalue_Mx512_My32-P010000-A0x400-A0y50-07-Aug-2010.mat'
%   %load 'MV9_Eigenvalue_Mx512_My32-P03000-A0x400-A0y50-05-Aug-2010.mat'
%   %load 'MV9_Eigenvalue_Mx512_My32-P05000-A0x400-A0y50-06-Aug-2010.mat'
% 
% elseif ( 0 )
%   
%   P0 = 10e3;
%   
%   load('BBN_z1.cm_ShermiX513_Npts513_282.84.mat');
%   hermiX = Expression1;
%   
%   load('BBN_z1.cm_ShermiX513_Npts513_35.355.mat');
%   hermiY = Expression1;
%     
% % load 'MV9_Eigenvalue_Mx128_My32-P03000-A0x200-A0y25-09-Aug-2010.mat'
% % load 'MV9_Eigenvalue_Mx128_My32-P05000-A0x200-A0y25-09-Aug-2010.mat'
%  load 'MV9_Eigenvalue_Mx128_My32-P010000-A0x200-A0y25-09-Aug-2010.mat'
% 
% % load 'MV9_Eigenvalue_Mx128_My32-P03000-A0x200-A0y50-09-Aug-2010.mat'
% % load 'MV9_Eigenvalue_Mx128_My32-P05000-A0x200-A0y50-09-Aug-2010.mat'  
% 
% 
% end

  %P0 = 20e3;
  P0 = 0.185e3;
  
  %load(sprintf('BBN_z0.cm_ShermiX513_Npts513_%g.mat',round(a0px/1e-6)));
  %load(sprintf('BBN_z0.cm_ShermiX2049_Npts513_%g.mat',round(a0px/1e-6)));
  %load(sprintf('3.2grid_BBN_z0.cm_ShermiX2049_Npts513_%g.mat',round(a0px/1e-6)));
  %load('UTA1.6mm_hermiX513_19.5.dat');
  
  %load('UTA_ShermiX513_Npts513_19.5.mat');%has size 0.6x0.6mm^2 grid
  load('UTA1.2mm_hermiX513_19.5.dat')
  hermiX = UTA1_2mm_hermiX513_19_5;
  %hermiX = Expression;
  %clear Expression
  
  %load(sprintf('BBN_z0.cm_ShermiX513_Npts513_%g.mat',round(a0py/1e-6)));
  %load(sprintf('BBN_z0.cm_ShermiX2049_Npts513_%g.mat',round(a0py/1e-6)));
  %load(sprintf('3.2grid_BBN_z0.cm_ShermiX2049_Npts513_%g.mat',round(a0py/1e-6)));
  %hermiY = Expression;
  hermiY = hermiX;  
  %clear Expression
  
  %load 'MV9_Eigenvalue_Mx512_My32-P010000-A0x400-A0y50-07-Aug-2010.mat'
  %load 'MV9_Eigenvalue_Mx512_My32-P03000-A0x400-A0y50-05-Aug-2010.mat'
  %load 'MV9_Eigenvalue_Mx512_My32-P05000-A0x400-A0y50-06-Aug-2010.mat'
  
  %this 12-May-2010 data set corresponds to n=2.14, x_eff = 2*d_eff points
  %load 'MV9_Eigenvalue_Mx2048_My08-deff-8.7-P05000-A0x800-A0y25-12-May-2011.mat'
  
  % this is actual BBN-value data set
  %load 'MV9_Eigenvalue_Mx2048_My08-deff-8.7-P05000-A0x800-A0y25-23-May-2011.mat'
  %load 'MV9_Eigenvalue_Mx2048_My08-deff-8.7-P010000-A0x800-A0y25-01-Jun-2011.mat'
  %load 'MV9_Eigenvalue_Mx2048_My08-deff-8.7-P015000-A0x800-A0y25-26-May-2011.mat'
  %load 'MV9_Eigenvalue_Mx2048_My08-deff-8.7-P020000-A0x800-A0y25-29-May-2011.mat'
    
  addpath('./Fourier-imaging/HarrisExp-071811/')
  load 'MV9_Eigenvalue_Mx32_My32-deff-8-P0185-A0x19.5-A0y19.5-18-Jul-2011.mat';

psastr = sprintf('P0=%g-kW-spot-%dx%d-Mx%d-My%d',P0/1e3,round(a0px/1e-6), round(a0py/1e-6),Mx,My);

l=flipud(l);
evec = fliplr(evec);
vecs = evec(1:Mx*My,1:nmodes) + 1i*evec(1+Mx*My:end,1:nmodes);

eigenmodes = zeros( 513, 513, nmodes );

%% visualize the eigenvector/mode
for idmode = 1:nmodes;
     fprintf('Working on Eigenmode %g\n',idmode);
     mode = zeros(size(hermiY,2),size(hermiX,2));
     idx = 0;

     mode = [hermiX(1:Mx,:).'*reshape(vecs(:,idmode),My,Mx).'*hermiY(1:My,:)].';

     %% all modes are orthonormal over the whole space, but in the
     %% limited grid size of 1.6x1.6mm^2 we make them orthogonal
     %% again.
     %%mode = mode/sqrt((sum(sum(abs(mode.^2)))*(1.6e-3/512).^2));
     mode_energy = (sum(sum(abs(mode.^2)))*dx*dy);
     fprintf('Energy of mode # %d = %g\n', idmode, mode_energy);
     mode = mode/sqrt(mode_energy);
     %%sanity check : modes are anyways orthonormal
     %%sum(sum(abs(mode.^2)))*(1.6e-3/512).^2
     %%sum(sum(abs(mode.^2)))*(3.2e-3/512).^2
     
     eigenmodes(:,:,idmode) = mode;

    if ( 1 )
     close all
     figure
     surf(x(1:2:end)/1e-3,y(1:2:end)/1e-3,abs(mode(1:2:end,1:2:end)).^2)
     view(0,90)
     shading interp
     title(sprintf('Mode# %d, [Eigval=%g] , %s\\mum^2',idmode-1,l(idmode),psastr))
     xlabel(['x axis (mm),', num2str(Mx),'HG modes']); 
     ylabel(['y axis (mm),', num2str(My),'HG modes']);
     set(gcf,'Color','White')
     axis square
     colorbar
     %pause(1)
     print(gcf,'-dpng',['BBN-z=0cm-Spatial-Eigenmode-',num2str(idmode-1),'-Eigenvalue=',num2str(l(idmode)),psastr,'.png']);
    end
end
save('-mat',['fast_xtalcenter_MV9_Eigenmodes_',psastr,'.mat'],'eigenmodes','l','evec')

% save('-mat',['BBN-z=1cm-MV9_Eigenmodes',num2str(nmodes),'_',psastr,'.mat'],'eigenmodes','l','evec','x')
save('-mat',['BBN-z=0cm-MV9_Eigenmodes',num2str(nmodes),'_',psastr,'.mat'],'eigenmodes','l','evec','x')
fprintf('%s\n',['Saved file : MV9_Eigenmodes',num2str(nmodes),'_',psastr,'.mat'])
