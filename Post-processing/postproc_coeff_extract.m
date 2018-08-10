% (C) 2011-2012, Muthiah Annamalai
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% 01/20/12
% --------
% Post-processing codes to extract phase information of the first 20
% eigenmodes via storing their complex coefficients in the Sqrt[2]
% calculation basis.
%
clear all
close all
addpath('./Scalefactor-data/')
N=20; %number of eigenmodes to save data from

%datafile = 'MV9_Eigenvalue_Mx32_My32-deff-8.7-P01150-A0x25-A0y25-10-Jun-2011';
%datafile = 'MV9_Eigenvalue_Mx32_My32-P02400-A0x100-A0y25-31-Mar-2010';
%datafile = 'MV9_Eigenvalue_Mx32_My32-P03100-A0x100-A0y50-20-Apr-2010';
%datafile = 'MV9_Eigenvalue_Mx32_My32-P05000-A0x100-A0y100-20-Apr-2010';
%datafile = 'MV9_Eigenvalue_Mx128_My32-P09000-A0x200-A0y100-15-Apr-2010';
%datafile = 'MV9_Eigenvalue_Mx128_My128-P016250-A0x200-A0y200-10-Apr-2010';
%datafile ='./MV9-EIGENVAL-DATA/MV9_Eigenvalue_Mx512_My32-P017250-A0x400-A0y100-06-Apr-2010.mat';
%datafile = './MV9-EIGENVAL-DATA/MV9_Eigenvalue_Mx2048_My08-P021500-A0x800-A0y50-25-Apr-2010.mat';

%only eigenvalue data no sign of eigenvectors!
%%datafile = 'MV9_Eigenvalue_Mx128_My128-P020000-A0x400-A0y100-27-Jan-2010';

Mx=32; My=32;
load([datafile,'.mat'])
l=flipud(l);
evec = fliplr(evec);
vecs = evec(1:Mx*My,1:N) + 1i*evec(1+Mx*My:end,1:N);
for idmode = 1:N
    hgmode = reshape(  vecs(:,idmode).^2, My, Mx ).';
    fprintf('sanity check %g ?= 1 \n',sum(sum(abs(hgmode))))
    save('-mat',sprintf('Eigenmode%d-AmpCoeff-Sqrt2-%s.mat',idmode-1,datafile))
end
