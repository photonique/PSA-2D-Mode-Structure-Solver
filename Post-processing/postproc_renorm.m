% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%% first 100 eigenvalues : proc 03/30/10
%% ordering of eigenvalue/vectors is from the largest to smallest.
%% This code tries to redo the images generated for CLEO-10
%% for eigenvalue spectrum and the  HG-mode distributions for the
%% most squeezed eigenvector.
clear all
close all
Mx = 32; My =32;
%load 'MV9_Eigenvalue_Mx32_My32-P01150-A0x25-A0y25-31-Mar-2010.mat'
load 'MV9_Eigenvalue_Mx32_My32-P05000-A0x100-A0y100-31-Mar-2010.mat'
%load 'MV9_Eigenvalue_Mx32_My32-P02400-A0x100-A0y25-31-Mar-2010.mat'
l=flipud(l);
plot(1./l,'-ob');
title('Eigenvalue spectrum 18.75kW pump with 400x100\mum^2 spot evaluated at 512x32 modes')
xlabel('Mode index')
ylabel('Eigenvalue')
set(gcf,'Color','White')

%% process smallest 15 eigenvectors
evec = fliplr(evec);
vecs = evec(1:Mx*My,1:15) + 1i*evec(1+Mx*My:end,1:15);

HGModes = zeros(Mx,My,15);
for idx = 1:1
 HGModes(:,:,idx)=reshape( abs( vecs(:,idx) ).^2, My, Mx )';
 figure
 surf(0:My-1,0:Mx-1,HGModes(:,:,idx));
 shading flat
 axis([0 My-1, 0 Mx-1, 0 max(max(HGModes(:,:,idx)))]);
 colorbar
 title(['Eigenvector=[',num2str(idx-1),',Eigenvalue=',num2str(l(idx)),'] { particular pump power, spot size with 32x32 modes}'])
 xlabel('Y-Mode index ')
 ylabel('X-Mode index ')
 view(0,90)
 set(gcf,'Color','White')
end
