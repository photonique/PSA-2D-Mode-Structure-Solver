%            
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
clear all
echo on
%% Mx = 256; My = 32;
P0 = 20e3; ap0x = 400e-6; ap0y = 100e-6;

%procfilename = {'proc-MV9_Mx256_My32-Wed-Jan-27-13-48-09-2010-eo.mat',
%	       'proc-MV9_Mx256_My32-Wed-Jan-27-13-48-09-2010-oe.mat',
%	       'proc-MV9_Mx256_My32-Wed-/home/muthu/Desktop/send-vJan-27-13-48-09-2010-ee.mat',
%	       'proc-MV9_Mx256_My32-Wed-Jan-27-13-48-09-2010-oo.mat'};

Mx = 512; My = 32;
procfilename = {'proc-MV9_Mx512_My32-Thu-Jan-28-08-44-41-2010-ee.mat',
		'proc-MV9_Mx512_My32-Thu-Jan-28-08-44-41-2010-eo.mat',
		'proc-MV9_Mx512_My32-Thu-Jan-28-08-44-41-2010-oe.mat',
		'proc-MV9_Mx512_My32-Thu-Jan-28-08-44-41-2010-oo.mat'};

for idx = 1:4
  procfile = procfilename{idx};
  load(procfile);
end

%% complex valued Greens function : QDs(1:Mx*My,:) + 1i*QDs(Mx*My+1:end,:);

%% inline this 
%% QDs = [blockset(QDse(:,1:ceil(Mx/2)*My),QDso(:,1:ceil(Mx/2)*My),Mx,My),blockset(QDse(:,(1+ceil(Mx/2)*My):end),QDso(:,(1+ceil(Mx/2)*My):end),Mx,My)];

QDs = zeros(Mx*My,2*Mx*My); QDse = zeros(Mx*My/2,Mx*My); QDso = QDse;

%% level-2 rec
QDse(1:2:end,1:2:end)=QDsee; QDso(1:2:end,1:2:end)=QDsoe;
QDse(2:2:end,2:2:end)=QDseo; QDso(2:2:end,2:2:end)=QDsoo;
clear QDsee QDseo QDsoe QDsoo

%% level-1 rec
QDs = [blockset(QDse(:,1:ceil(Mx/2)*My),QDso(:,1:ceil(Mx/2)*My),Mx,My),blockset(QDse(:,(1+ceil(Mx/2)*My):end),QDso(:,(1+ceil(Mx/2)*My):end),Mx,My)];
clear QDse QDso

tic;
%% calculate the eigenmodes 'v', and squeezing spectrum 'l'.
QDs = [real(QDs); imag(QDs)];
G = full(QDs*(QDs.'));
opts.isreal = 1; opts.issym = 1;
%% [evec,l] = eig(G);
[evec,l] = eigs( G,100,'sm',opts );
l = diag(l);
time_lapse = toc

flipud(l)
