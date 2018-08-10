% (C) 2011 Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
%% 10/08/10 : Some minor changes to grid limits (x_top,x_bot) to make
%% things smoother
%% 10/07/10 : Use grid for NU-SLM experiment over 190x190 um^2 square grid
%% first eigenmodes : proc 03/30/10
%% ordering of eigenvalue/vectors is from the largest to smallest.
%% visualize the eigenmodes properly in spatial domain
%% 04/1/10 : back to mathematica generated polynomials.
%% 05/03/10: Change for plotting on the smaller scale from +/-1.2mm to
%% +/-0.3mm
function spatialmodevis(Mx,My,P0,a0px,a0py,vecs,l,dosave)
if ( nargin < 8 )
    dosave =  false;
end
%% Standard HG-mode calculations over 1.6x1.6 mm^2 square grid at 501x501 points
% Xpts = 1024; Ypts = Xpts; x = linspace(-1.6e-3,+1.6e-3,Xpts); y = x;
% dx = 3.2e-3/(Xpts-1); dy = 3.2e-3/(Ypts-1);
% fprintf('Spatial modevis: 1024x1024 pt, 1.6x1.6mm^2 data set\n')
% hermiX = load(sprintf('UTA3.2mm_hermiX513_%d.dat',round(a0px/1e-6)));
% hermiY = load(sprintf('UTA3.2mm_hermiX513_%d.dat',round(a0py/1e-6)));
Xpts = 501; 513; 1024; Ypts = Xpts; x = linspace(-1.6e-3,+1.6e-3,Xpts); y = x;
dx = 3.2e-3/(Xpts-1); dy = 3.2e-3/(Ypts-1);
%fprintf('Spatial modevis: 1024x1024 pt, 1.6x1.6mm^2 data set\n')
%hermiX = load(sprintf('UTA1.6mm_hermiX513_%d.dat',round(a0px/1e-6)));
%hermiY = load(sprintf('UTA1.6mm_hermiX513_%d.dat',round(a0py/1e-6)));

%load('UTA1.6mm_hermiX513_19.5.dat');
%hermiX = UTA1_6mm_hermiX513_19_5;
load('UTA1.6mm_hermiX513_25.dat');
hermiX = UTA1_6mm_hermiX513_25;
hermiY = hermiX;
%load('ShermiX513_27.dat');
%hermiX = ShermiX513_27;
%load('ShermiX513_25.7.dat')
%hermiY = ShermiX513_25_7;
%
%load('UTA1.6mm_hermiX513_19.5.dat');
%hermiX = UTA1_6mm_hermiX513_19_5;
%hermiY = hermiX;
%% Standard HG-mode calculations over 3.2x3.2 mm^2 square grid at 501x501 points
% Xpts = 501; Ypts = Xpts; x = linspace(-1.6e-3,+1.6e-3,Xpts); y = x;
%fprintf('Spatial modevis: 501x501 pt, 3.2x3.2mm^2 data set\n')
%hermiX = load(sprintf('hermiX2048_%d.dat',round(a0px/1e-6)));
%hermiY = load(sprintf('hermiX2048_%d.dat',round(a0py/1e-6)));
%hermiX = load(sprintf('hermiX513_%d.dat',round(a0px/1e-6)));
%hermiY = load(sprintf('hermiX513_%d.dat',round(a0py/1e-6)));

%% Mode visualization for BBN data sets from Aug 2010
%fprintf('Spatial modevis: BBN data set\n')
%Xpts = 513; Ypts = 513; x = linspace(-1.6e-3/2, +1.6e-3/2, Xpts); y = x;
% load(sprintf('BBN_z0.cm_ShermiX513_Npts513_%d.mat',round(a0px/1e-6)));
% hermiX = Expression;
% load(sprintf('BBN_z0.cm_ShermiX513_Npts513_%d.mat',round(a0py/1e-6)));
% hermiY = Expression;

%% NU-SLM experiment over 190x190 um^2 square grid with 512 x 512 points.
%Xpts = 512; Ypts = 512; x = linspace(-371.0938e-6/2, +371.0938e-6/2, Xpts); y = x;
%fprintf('Spatial modevis: NU SLM 512x512 point, 371x371um^2 grid data set\n')
%hermiX = load(sprintf('NU371um_hermiX36_%d.dat',round(a0px/1e-6)));
%hermiY = load(sprintf('NU371um_hermiX36_%d.dat',round(a0py/1e-6)));
%hermiX = load('NU371um_hermiX36_36.9.dat');
%hermiY = load('NU371um_hermiX36_36.9.dat');
%hermiX = load(sprintf('NU190um_hermiX36_%d.dat',round(a0px/1e-6)));
%hermiY = load(sprintf('NU190um_hermiX36_%d.dat',round(a0py/1e-6)));

a0sx = sqrt(2)*a0px; a0sy = sqrt(2)*a0py;
%% vecs = evec(1:Mx*My,1:15) + 1i*evec(1+Mx*My:end,1:15);
%% visualize the eigenvector/mode 0
mode = zeros(size(hermiY,2),size(hermiX,2));
EM = [];HGModes = [];
if ( dosave )
    mkdir('HG-mode-data')
    mkdir('HG-mode-profile')
    mkdir('spatial-eigenmode-data')
    mkdir('Spatial-Eigenmode-images')
end
%y_top = find( x > 0.3e-3,1 ); y_bot = find( x > -0.3e-3,1 );
%x_top = find( x > 1.2e-3,1 ); x_bot = find( x > -1.2e-3,1 );
%x = x(x_bot:x_top)/1e-3; y = y(y_bot:y_top)/1e-3;
%y_top = 513; y_bot = 1; %BBN data set only
%y_top = 512; y_bot = 1; %generic processing 
%y_top = 1024; y_bot = 1; %Hresolution data set
%x_top = y_top; x_bot = y_bot;
for idmode = 1:size(vecs,2);
    fprintf('Working on Eigenmode %g\n',idmode);
    hgmode = reshape( abs( vecs(:,idmode) ).^2, My, Mx )';
    HGModes(1:Mx,1:My,idmode) = hgmode;
    size(hermiX(1:Mx,:))
    size(vecs(:,idmode))
    size(hermiY(1:My,:))
    mode = hermiX(1:Mx,:)'*reshape(vecs(:,idmode),My,Mx)'*hermiY(1:My,:);    
    figure
    %mode = mode(x_bot:x_top,y_bot:y_top);
    mode = mode.';
    fprintf('Energy in mode # %d = %g \n',idmode-1,sum(sum(abs(mode).^2))*dx*dy);    
    modeint = abs(mode).^2;
    %size(mode); size(modeint)
    % EM(idmode,1:Xpts,1:Ypts) = modeint;
    % surf(x/1e-3,y/1e-3,modeint)
    surf(x/1e-3,y/1e-3,modeint)
    shading interp
    title(sprintf('Mode# %d, [Eigval=%g] , P0=%gk/spot:%dx%d\\mum^2',idmode-1,l(idmode),P0/1e3,round(a0px/1e-6),round(a0py/1e-6)))
    xlabel(sprintf('x axis (mm), %d HG modes',Mx)); ylabel(sprintf('y axis (mm), %d HG modes',My));    
    set(gcf,'Color','White')
    % axis([[-1.6e-3,+1.6e-3, -1.6e-3,1.6e-3]/1e-3,0,max(max(modeint))])
    %axis([[-1.2e-3,+1.2e-3, -0.3e-3,0.3e-3]/1e-3,0,max(max(modeint))])
    axis equal
    colorbar
    view(0,90)
    if ( dosave )
      print(gcf,'-dpng',['./Spatial-Eigenmode-images/Spatial-Eigenmode-',num2str(idmode-1),'-Eigenvalue=',num2str(l(idmode)),sprintf('P0=%gkW-spot-%dx%d.png',P0/1e3,round(a0px/1e-6),round(a0py/1e-6))]);
      save('-v6',['./spatial-eigenmode-data/Spatial-Eigenmode-',num2str(idmode-1),'-Eigenvalue=',num2str(l(idmode)),sprintf('P0=%gkW-spot-%dx%d.mat',P0/1e3,round(a0px/1e-6),round(a0py/1e-6))],'mode','modeint');
      save('-v6',['./HG-mode-data/Mode-',num2str(idmode-1),'-Eigenvalue=',num2str(l(idmode)),sprintf('P0=%gkW-spot-%dx%d.mat',P0/1e3,round(a0px/1e-6),round(a0py/1e-6))],'hgmode');    
    end
    figure
    surf([0:(Mx-1)]',[0:(My-1)],hgmode.')
    shading flat
    title(sprintf('|Amn|^2 Mode# %d, [Eigval=%g] , P0=%gk/spot:%dx%d\\mum^2',idmode-1,l(idmode),P0/1e3,round(a0px/1e-6),round(a0py/1e-6)))
    xlabel(sprintf('x modes , %d HG modes',Mx)); ylabel(sprintf('y modes, %d HG modes',My));
    set(gcf,'Color','White')
    axis([0 Mx-1,0,My-1,0,1])
    axis equal
    colorbar
    view(0,90)    
    if ( dosave )        
       print(gcf,'-dpng',['./HG-mode-profile/HG-Eigenmode-',num2str(idmode-1),'-Eigenvalue=',num2str(l(idmode)),sprintf('P0=%gkW-spot-%dx%d.png',P0/1e3,round(a0px/1e-6),round(a0py/1e-6))]);
    end
end

