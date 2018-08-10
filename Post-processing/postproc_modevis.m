% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%% 
%% 08/04/10 : Looks like the mode absolute value itself is
%% hermitian-conjugated accidentally while being transposed;
%% this may not be a problem for visualizing the eigenmodes at z=0, but
%% maybe later. So we resort to using just simple transpose.
%%
%% 04/11/10 04/12/10 : we try visualized eigenmodes of 200x200 spot size 
%%                     with pump power 16.25kW for 128x128 polynomials.
%% 
%% first eigenmodes : proc 03/30/10
%% ordering of eigenvalue/vectors is from the largest to smallest.
%% visualize the eigenmodes properly in spatial domain
%% 04/1/10 : back to mathematica generated polynomials.
%% 

clear all
close all
more off
%% Xpts = 1500; Ypts = 150; x = linspace(-10e-3,+10e-3,Xpts); y = linspace(-1e-3,+1e-3,Ypts);
Xpts = 501; Ypts = Xpts; x = linspace(-1.6e-3,+1.6e-3,Xpts); y = x;
%Mx = 512; My = 32; a0px = 400e-6; a0py = 100e-6;
Mx = 128; My = 128; a0px = 200e-6; a0py = 200e-6;
a0sx = sqrt(2)*a0px; a0sy = sqrt(2)*a0py;
%hermiX = load('hermiX513_400.dat');
%hermiY = load('hermiX513_100.dat');
hermiX = load('hermiX513_200.dat');
hermiY = hermiX;
%load 'MV9_Eigenvalue_Mx512_My32-P018750-A0x400-A0y100-26-Mar-2010.mat'
%load 'MV9_Eigenvalue_Mx512_My32-P017250-A0x400-A0y100-06-Apr-2010.mat'
load 'MV9_Eigenvalue_Mx128_My128-P016250-A0x200-A0y200-10-Apr-2010.mat'
l=flipud(l);
evec = fliplr(evec);
vecs = evec(1:Mx*My,1:15) + 1i*evec(1+Mx*My:end,1:15);

%% form the HermiX,Y from normhermitepolynomial.
%%hermiX = zeros(Mx,Xpts); hermiY = zeros(My,Ypts);
% for mx = 1:Mx
%     hermiX(mx,:) = normhermitepoly(mx-1,x/a0sx).*exp(-0.5*(x/a0sx).^2)/a0sx;
% end
% for my = 1:My
%     hermiY(my,:) = normhermitepoly(my-1,y/a0sy).*exp(-0.5*(y/a0sy).^2)/a0sy;
% end

%% visualize the eigenvector/mode 0
mode = zeros(size(hermiY,2),size(hermiX,2));
size(mode)
%figure
for idmode = 1:15;
    fprintf('Working on Eigenmode %g\n',idmode);
    mode = hermiX(1:Mx,:).'*reshape(vecs(:,idmode),My,Mx).'*hermiY(1:My,:);
    figure
    %subplot(5,3,idmode)
    surf(x/1e-3,y/1e-3,abs(mode').^2)
    shading interp
    %title(sprintf('Mode# %d, [Eigval=%g] , P0=18.75k/spot:400x100\\mum^2',idmode-1,l(idmode)))
    %title(sprintf('Mode# %d, [Eigval=%g] , P0=17.25k/spot:400x100\\mum^2',idmode-1,l(idmode)))
    title(sprintf('Mode# %d, [Eigval=%g] , P0=16.25k/spot:200x200\\mum^2',idmode-1,l(idmode)))
    %xlabel('x axis (mm), 512 HG modes'); ylabel('y axis (mm),  32 HG modes');
    xlabel('x axis (mm), 128 HG modes'); ylabel('y axis (mm),  128 HG modes');
    set(gcf,'Color','White')
    %axis([[-1.6e-3,+1.6e-3, -0.5e-3,0.5e-3]/1e-3,0,max(max(abs(mode).^2))])
    axis([[-0.5e-3,+0.5e-3, -0.5e-3,0.5e-3]/1e-3,0,max(max(abs(mode).^2))])
    axis square
    %axis tight
    colorbar
    view(0,90)
    print(gcf,'-dpng',['Spatial-Eigenmode-',num2str(idmode-1),'-Eigenvalue=',num2str(l(idmode)),'-16.25kW-spot200x200.png']);
end
%     idx = 0;
%     for my = 1:32
%         for mx = 1:512
%             idx = idx + 1;
%             mode = mode + vecs(idx,idmode)*(hermiY512_100(my,:).')*(hermiX512_400(mx,:));
%         end
%     end
