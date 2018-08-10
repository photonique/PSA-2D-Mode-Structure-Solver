% (C) 2011 Muthiah Annamalai
%
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%% Driver script for overlap integral calculation
%%load 'Spatial-Eigenmode-0-Eigenvalue=0.063283P0=9.25kW-spot-200x100.mat'
%%a0px = 200e-6/sqrt(10); a0py = 100e-6/sqrt(10);
%  
%  
%  05/18/10 : Calculate the overlap integral for higher order elliptic modes
%  with a specified spatial eigenmode. By default we only seek match to 0-0
%  mode still.
%  
clear all
close all

mode_m = 0; mode_n = 0;
%% load './Scalefactor-data/Spatial-Eigenmode-0-Eigenvalue=0.06557P0=17.25kW-spot-400x100.mat'
%% best spot size was 84.8485 x 48.4848 with overlap 0.992557
%% best spot size was 84.8485 x 48.4848 with overlap 0.992557
%% load './Scalefactor-data/Spatial-Eigenmode-0-Eigenvalue=0.065634P0=2.4kW-spot-100x25.mat'
%% best spot size was 51.5152 x 31.8182 with overlap 0.984098
%% best spot size was 51.5152 x 31.4394 with overlap 0.984106
%% load './Scalefactor-data/Spatial-Eigenmode-0-Eigenvalue=0.065785P0=9kW-spot-200x100.mat'
%% best spot size was 62.2449 x 49.4898 with overlap 0.991229
%% best spot size was 63.6364 x 48.4848 with overlap 0.991479
%% load './Scalefactor-data/Spatial-Eigenmode-0-Eigenvalue=0.065802P0=3.1kW-spot-100x50.mat'
%% best spot size was 51.5306 x 40.051 with overlap 0.984414
%% best spot size was 50 x 39.3939 with overlap 0.984534

addpath('Scalefactor-data/')
% load 'Spatial-Eigenmode-0-Eigenvalue=0.0011902P0=16.25kW-spot-100x50.mat'
% best spot size was 43.9394 x 37.1212 with overlap 0.90367
% load 'Spatial-Eigenmode-0-Eigenvalue=0.00030895P0=16.25kW-spot-100x25.mat'
% best spot size was 42.4242 x 31.8182 with overlap 0.864071
% load 'Spatial-Eigenmode-0-Eigenvalue=0.18125P0=1.01562kW-spot-100x25.mat'
% best spot size was 54.5455 x 31.8182 with overlap 0.991957
% load 'Spatial-Eigenmode-0-Eigenvalue=0.11394P0=2.03125kW-spot-100x50.mat'
% best spot size was 51.5152 x 40.1515 with overlap 0.989806
% best spot size was 52.5253 x 40.404 with overlap 0.98969

%load 'Spatial-Eigenmode-0-Eigenvalue=0.071185P0=16.25kW-spot-400x100.mat'
%A0PX = 400e-6; A0PY = 100e-6;
% overlap-400x100-p016-25kW-v1.png
% best spot size was 100 x 47.9798 with overlap 0.976238
% best spot size was 81.6327 x 48.9796 with overlap 0.992912
% load 'Spatial-Eigenmode-0-Eigenvalue=0.092294P0=16.25kW-spot-800x50.mat'
% best spot size was 113.131 x 38.3838 with overlap 0.99005
% load 'Spatial-Eigenmode-8-Eigenvalue=0.11243P0=21.5kW-spot-800x50.mat'
% A0PX = 800e-6; A0PY = 50e-6;

%load 'Spatial-Eigenmode-0-Eigenvalue=0.063622P0=21.5kW-spot-800x50.mat'
% best spot size was 114.286 x 38.7755 with overlap 0.991804

% find best overlap for spatial eigenmode 4 with HG40 mode.
%load 'Spatial-Eigenmode-4-Eigenvalue=0.082877P0=21.5kW-spot-800x50.mat'
addpath('./Fourier-imaging')
%load 'Spatial-Eigenmode-0-Eigenvalue=0.069943P0=6.5kW-spot-800x50.mat'

load 'Spatial-Eigenmode-5-Eigenvalue=0.093669P0=6.5kW-spot-800x50.mat'
% eigenmode 5 saving only the last side-lobe
mode_5side = mode;
mode_5side(:,1:600) = 0;

% recenter this mode at 512,512
% from its present (512,622) which is 0.38162mm from (0,0) in grid.
rc_5side = [mode_5side(:,122:end), zeros(1024,121)];
mode = rc_5side;

mode_m = 0; mode_n = 0;
A0PX = 800e-6; A0PY = 50e-6;

tic
ovrlap = [];
SF = linspace(0.05,0.5,25);
for idx = 1:length(SF)
    for idy = 1:length(SF)
      sfx = SF(idx); sfy = SF(idy);
      a0px = A0PX*sfx/5; a0py = A0PY*sfy;
      A0px(idx) = a0px; A0py(idy) = a0py;
      ovrlap(idx,idy) = spotsize_scalefactor(a0px,a0py,mode,mode_m,mode_n);
    end
end
toc
surf(A0py/1e-6,A0px/1e-6,ovrlap)
shading flat
colorbar
axis([0, A0py(end)/1e-6, 0, A0px(end)/1e-6, 0, 1])
view(0,90)
[idx,idy]= find(ovrlap >= max(max(ovrlap)) );
fprintf('best spot size was %g x %g with overlap %g\n',A0px(idx)/1e-6,A0py(idy)/1e-6,ovrlap(idx,idy));
title(sprintf('best spot size was %g x %g \\mum^2 with overlap %g',A0px(idx)/1e-6,A0py(idy)/1e-6,ovrlap(idx,idy)))
xlabel('Y spot size sig \mum')
ylabel('X spot size sig \mum')
zlabel('overlap integral with 0-eigen mode')
set(gcf,'Color','White')

% values of optimum data for the elliptic pump
% from studies in 04/20/10.
%    pump X    optim X    pump Y   optim Y
%   100.0000   51.5152   25.0000   31.8182
%   100.0000   51.5152   50.0000   40.0510
%   200.0000   62.2449  100.0000   49.4898
%   400.0000   84.8485  100.0000   48.4848
% 

% inlcuded data from round-spot size also.
% figure
% pumpspotX  = [25, 100, 100, 100, 200,200,  400, 800 ];
% optimspotX = [32.1572, 49.5495, 51.5152, 51.5152, 62.2449, 62.0621, 84.8485, 113.131 ];
% pumpspotY  = [25, 25, 50, 50, 100,  100, 100, 200  ];
% optimspotY = [32.1572, 31.8182, 38.3838, 40.051, 49.4898, 48.4848 , 49.5495, 62.0621];
% 
% %% normalized intensity 0.40625w/um^2.
% ni_pumpX = [25, 100, 100, 100];
% ni_optimX = [32.78280	54.54550	52.52530	50.05010];
% ni_pumpY = [25	25	50	100];
% ni_optimY = [32.78280	31.18182	40.40400	50.05010];
% 
% %% power set at 16.25kW
% p200_pumpX = [25	100	100	100];
% p200_optimX = [32.78280	43.93900	42.42420	45.54550];
% p200_pumpY = [25	25	50	100];
% p200_optimY = [32.78280	31.18182	37.12120	45.54550];
% 
% loglog(pumpspotX/25,optimspotX/25,'-ob',pumpspotY/25,optimspotY/25,'-or',...
%    ni_pumpX/25,ni_optimX/25,'-og',ni_pumpY/25,ni_optimY/25,'-ok',...
%    p200_pumpX/25,p200_optimX/25,'-om',p200_pumpY/25,p200_optimY/25,'-oy')
% legend({'X-spot size @ gain 15','Y-spot size @ gain 15','X-spot @ norm intensity','Y-spot @ norm intensity',...
%     'X-spot size @ P0=16.25kW','Y-spot size @ P0=16.25kW'})
% xlabel('pump spot / 25\mum, in log10-scale')
% ylabel('signal optimum spot/ 25\mum, in log10-scale')
% title('Optimum value of the signal overlap spot size for given Elliptic and Round pump spot')
% set(gcf,'Color','White')
