% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%% first 100 eigenvalues : proc 03/30/10
%% ordering of eigenvalue/vectors is from the largest to smallest.
clear all
close all
Mx = 512; My =32;
%load 'MV9_Eigenvalue_Mx512_My32-P018750-A0x400-A0y100-26-Mar-2010.mat'
%load 'MV9_Eigenvalue_Mx512_My32-P017250-A0x400-A0y100-06-Apr-2010.mat'
load 'MV9_Eigenvalue_Mx128_My128-P016250-A0x200-A0y200-10-Apr-2010.mat'
l=flipud(l);
if ( 1 )
semilogy(l,'-ob');
%title('Eigenvalue spectrum 18.75kW pump with 400x100\mum^2 spot evaluated at 512x32 modes')
%title('Eigenvalue spectrum 17.25kW pump with 400x100\mum^2 spot evaluated at 512x32 modes')
title('Eigenvalue spectrum 16.25kW pump with 200x200\mum^2 spot evaluated at 128x128 modes')
xlabel('Mode index')
ylabel('Eigenvalue')
set(gcf,'Color','White')
end
%% process smallest 15 eigenvectors
evec = fliplr(evec);
vecs = evec(1:Mx*My,1:15) + 1i*evec(1+Mx*My:end,1:15);

HGModes = zeros(Mx,My,15);
figure
for idx = 1:15
 HGModes(:,:,idx)=reshape( abs( vecs(:,idx) ).^2, My, Mx )';
 Mode = HGModes(:,:,idx);
 save('-v6',['Mode',num2str(idx-1),'-16.25kW-200x200umsq-041110.mat'],'Mode');
 figure
 %subplot(5,3,idx)
 surf(0:My-1,0:Mx-1,HGModes(:,:,idx)),shading flat;
 axis([0 My-1, 0 Mx-1, 0 max(max(HGModes(:,:,idx)))]);view(0,90);
 colorbar
 title(['Eigenvector=[',num2str(idx-1),',Eigenvalue=',num2str(l(idx)),'{ 16.25kW pump 200x200\mum^2 spot 128x128 modes}']);
 xlabel('Y-Mode index ')
 ylabel('X-Mode index ')
 set(gcf,'Color','White')
 print(gcf,'-dpng',['Interp-FlatTop-Eigenmode-',num2str(idx-1),'-Eigenvalue=',num2str(l(idx)),'-16.25kW-spot200x200.png']);
end
