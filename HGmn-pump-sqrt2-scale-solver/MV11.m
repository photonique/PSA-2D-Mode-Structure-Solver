% 
% (C) 2009, 2010 Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%     Nonlinear Optics & Nanophotonics Lab, UT-Arlington.
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% 07/12/10 : Use the MV9 template to write code for higher-order pump
%            mode eigenvalue solver, where the signal spot = sqrt(2) pump
%            scale is used.
%            
%            Essentially the overlap integral factorizes as z=0 part and
%            a z-dependent complex phase that is  easily calculated.
%            
%            Still the ways of calculating the overlap factors are subject
%            to numerical overflow. Resorted to older way.
% 
% 07/16/10 : Finishing up the whole solver. Breathe in life.
% 07/18/10 : Move scale factor into the Kappa.
% 07/18/10 : Begin to use the Mathematica generated data files.
%            Also importantly I think code is generalized for even pump modes only.
%            i.e : q%2=m%2 = 0
%            Setup the studies for the giant case 200x200 um^2 @ 128x128
%            modes for the pump order 4,4.
%            
clear all
close all
more off
tic
%% Number of signal-idler modes.
%%Mx = 16; My = 16;
%% Mx = 16; My = 16;
%% Mx = 4; My = 4;
%Mx = 32; My = 32;
%Mx = 16; My = 16;
%Mx = 6; My = 4;
%Mx = 4; My = 6;
%Mx = 2; My = 2;
%Mx = 32; My = 32;
%%Mx = 64; My = 32;
%%Mx = 1024; My = 32;
%Mx = 512; My = 32;
%Mx =256; My = 64;
Mx = 128; My = 128;
%Mx = 2048; My = 8;
%%Mx = 128; My = 32;
dosave = 1;
if ( mod(Mx,2) + mod(My,2) ~= 0 )
  error(' Mx, My need to be even simultaneously ')
end

%% using MKS units
c=2.99792458e8; eps0=8.85e-12;

%% NU-experiment parameters
%% Lz : 20mm crystal, signal: 1560nm, pump: 780nm, deff : 8.7e-12, ns,p: 1.78
lambda_s = 1560e-9; lambda_p = lambda_s/2; ks = 2*pi/lambda_s;  kp = ks*2;
Lz =  20e-3; ns = 1.78; np = ns; deff = 8.7e-12; omega_s = ks*c;
dk = 4*pi*(np-ns)/lambda_s;
theta_p = -pi/2;


% pump order : clearly we require 'even' order pump
% since we have the rule, m + l + n is always even for non-zero
% coupling integrals
% q = 2; m = 0;
q = 4; m = 4;

% we only handle the even order pump mode with this solver
assert( mod(q,2) + mod(m,2) == 0 )

%% simulation parameters
Lpts = 300; dz = Lz/(Lpts);

%% pump parameters: 20kW pump case with 1/e waist size of pump
%%P0 = 17.25e3; ap0x = 400e-6; ap0y = 100e-6;
%%P0 = 16.5e3; ap0x = 200e-6; ap0y = 200e-6;
%%P0 = 16.00e3; ap0x = 800e-6; ap0y = 50e-6;
%%P0 = 10.00e3; ap0x = 200e-6; ap0y = 100e-6;
%%P0 = 8.75e3; ap0x = 200e-6; ap0y = 100e-6;
%%P0 = 9.00e3; ap0x = 200e-6; ap0y = 100e-6;
%%P0 = 19.50e3; ap0x = 800e-6; ap0y = 50e-6;
%%P0 = 16e3; ap0x = 800e-6; ap0y = 50e-6;

%% Calibration with data for CLEO-10, and PQE-10.
%% Powers adjusted for uniform gain of 15.
%%P0 = 1.15e3; ap0x = 25e-6; ap0y = 25e-6;
%%P0 = 2.4e3; ap0x = 100e-6; ap0y = 25e-6;
%%P0 = 3.1e3; ap0x = 100e-6; ap0y = 50e-6;
%%P0 = 5e3; ap0x = 100e-6; ap0y = 100e-6;
%%
P0 = 16.25e3; ap0x = 200e-6; ap0y = 200e-6;
%%P0 = 21.5e3; ap0x = 800e-6; ap0y = 50e-6;
%%P0 = 17.25e3; ap0x = 400e-6; ap0y = 100e-6;

%% Calibrated for 16.25kW where 200x200 has gain 15.
% P0 = 16.25e3; ap0x = 25e-6; ap0y = 25e-6;
%%P0 = 16.25e3; ap0x = 100e-6; ap0y = 25e-6;
%%P0 = 16.25e3; ap0x = 100e-6; ap0y = 50e-6;
%%P0 = 16.25e3; ap0x = 100e-6; ap0y = 100e-6;
%%P0 = 16.25e3; ap0x = 800e-6; ap0y = 50e-6;

%% Calibrated for same intensity as 16.25kW where 200x200 has gain 15.
%% which is 0.40625 W/um^2.
%%P0 = 253.90625; ap0x = 25e-6; ap0y = 25e-6;
% P0 = 1015.62500; ap0x = 100e-6; ap0y = 25e-6;
% P0 = 2031.25; ap0x = 100e-6; ap0y = 50e-6;
% P0 = 4062.5; ap0x = 100e-6; ap0y = 100e-6;
% P0 = 21.5e3; ap0x = 800e-6; ap0y = 50e-6;
% P0 = 16.25e3; ap0x = 400e-6; ap0y = 100e-6;


fprintf('*************************************************************\n')
if ( dosave )
   fprintf('SAVE enabled\n')
else
   fprintf('Save DISABLED!\n')
end
fprintf('Running case : P0 = %g [kW], spot : %d x %d at modes %d x %d \n',P0/1e3, ap0x/1e-6,ap0y/1e-6,Mx,My);
fprintf('Pump mode q=%d, m=%d\n',q,m);
fprintf('*************************************************************\n')

as0x = sqrt(2)*ap0x; as0y = sqrt(2)*ap0y;

%% Rayleigh length for X, Y modes
Zrx = 2*pi*np*ap0x^2/lambda_p; Zry = 2*pi*np*ap0y^2/lambda_p;
Kappa = omega_s*deff/sqrt(2*eps0*ns^2*np*c^3)*sqrt(P0/(ap0x*ap0y));

%% fast Gouy phase calculation including a higher order pump
evenMx = 0.5 + (0:2:2*Mx) - q;
E0xe = kron(hankel(evenMx(1:Mx/2),evenMx(Mx/2:Mx-1)),toeplitz(mod(1:My,2)));
E0xo = kron(2 + hankel(evenMx(1:Mx/2),evenMx(Mx/2:Mx-1)),toeplitz(mod(1:My,2)));
oddMy = (0.5 + (0:(2*My-1)) - m).*mod([1:My*2],2);
E0ye = kron( ones(Mx/2,Mx/2), hankel( oddMy(1:My), oddMy(My:2*My-1)));
E0yo = E0ye;

E0xee = E0xe(1:2:end,1:2:end); E0xeo = E0xe(2:2:end,2:2:end);
E0xoe = E0xo(1:2:end,1:2:end); E0xoo = E0xo(2:2:end,2:2:end);

E0yee = E0ye(1:2:end,1:2:end); E0yeo = E0ye(2:2:end,2:2:end);
E0yoe = E0yo(1:2:end,1:2:end); E0yoo = E0yo(2:2:end,2:2:end);

clear E0ye E0yo E0xo E0xe evenMx oddMy

%% 
%% Indexing Conventions for 4D ((p,l),(r,n)) -> 2D (i,j) contraction 
%% notation. In new index conventions, i is row index, and j is col index,
%% which represent how the j-th mode couples to i-th mode.
%% 
%% Row : i == (p,l) mode  p : row index, l : col index
%%       p = mod( i, My ); l = floor( i / My );
%% Col : j == (r,n) modsparsee  transpose of ( r : row index, n : col index )
%%       r = mod( j, My ); n = floor( j / My );
%% keep everything row-vectors
%% 
pP = mod( [ 0:(Mx*My - 1)], My ); lL = floor( [ 0:(Mx*My - 1)] / My );
nN = lL; rR = pP;

LL = repmat( lL, [Mx*My,1] ); NN = LL';
%lambda_mln = calc_higherorder_overlap( q, Mx); %ap0x
lambda_mln = load_higherorder_overlap( q, Mx );
DLN = kron(lambda_mln,ones(My));
[DLNe,DLNo] = blockget(  DLN, Mx, My );
DLNee = DLNe(1:2:end,1:2:end); DLNeo = DLNe(2:2:end,2:2:end);
DLNoe = DLNo(1:2:end,1:2:end); DLNoo = DLNo(2:2:end,2:2:end);
clear NN LL DLN DLNe DLNo

PP = repmat( pP', [1, Mx*My] ); RR = PP';
% lambda_qpr = calc_higherorder_overlap( m, My); %ap0y
lambda_qpr = load_higherorder_overlap( m, My );
DPR = kron(ones(Mx),lambda_qpr);
[DPRe,DPRo] = blockget(  DPR, Mx, My );
DPRee = DPRe(1:2:end,1:2:end); DPReo = DPRe(2:2:end,2:2:end);
DPRoe = DPRo(1:2:end,1:2:end); DPRoo = DPRo(2:2:end,2:2:end);
clear RR PP DPR DPRe DPRo
clear pP lL nN rR DLNe DLNo

%% coupling integral is loop constant and we pull it out.
C0zee = full(Kappa*DPRee.*DLNee);
C0zoe = full(Kappa*DPRoe.*DLNoe);
C0zeo = full(Kappa*DPReo.*DLNeo);
C0zoo = full(Kappa*DPRoo.*DLNoo);

clear DPRoe DPRoo DPRee DPReo
clear DLNoo DLNoe DLNeo DLNee

%% calculate the Greens function
%% Excite the mn mode in both quadratures so that the
%% output contains the real-valued Greens function automatically.
%% equivalent to eye(2*Mx*My)

%% introducing even, odd blocks.
SZ2 = [Mx*My/4,Mx*My/4];

QDsee = [eye(SZ2),1i*eye(SZ2)]; 
QDseo = QDsee; QDsoe = QDsee; QDsoo = QDsee;

%% RK-4 method following Gilbert Strang
QDsn1ee = zeros(SZ2); QDsn2ee = QDsn1ee; QDsn3ee = QDsn1ee; QDsn4ee = QDsn1ee;
QDsn1eo = zeros(SZ2); QDsn2eo = QDsn1ee; QDsn3eo = QDsn1ee; QDsn4eo = QDsn1ee;
QDsn1oo = zeros(SZ2); QDsn2oo = QDsn1oo; QDsn3oo = QDsn1oo; QDsn4oo = QDsn1oo;
QDsn1oe = zeros(SZ2); QDsn2oe = QDsn1oo; QDsn3oe = QDsn1oo; QDsn4oe = QDsn1oo;

z = -Lz/2; zpts = 0;

while ( z < Lz/2 & zpts < Lpts )
  
  %% calculate Dsn1
  phizx = atan(z/Zrx); phizy = atan(z/Zry);
  scale_z = sqrt(sqrt(1+(z/Zrx)^2)*sqrt(1+(z/Zry)^2));
  Czee = 1i*C0zee.*exp(+1i*(theta_p+dk*z+E0xee*phizx + E0yee*phizy))./scale_z;
  Czeo = 1i*C0zeo.*exp(+1i*(theta_p+dk*z+E0xeo*phizx + E0yeo*phizy))./scale_z;
  Czoe = 1i*C0zoe.*exp(+1i*(theta_p+dk*z+E0xoe*phizx + E0yoe*phizy))./scale_z;
  Czoo = 1i*C0zoo.*exp(+1i*(theta_p+dk*z+E0xoo*phizx + E0yoo*phizy))./scale_z;
  
  QDsn1ee = Czee*conj(QDsee)./2.0;
  QDsn1eo = Czeo*conj(QDseo)./2.0;
  QDsn1oe = Czoe*conj(QDsoe)./2.0;
  QDsn1oo = Czoo*conj(QDsoo)./2.0;
  
  %% calculate Dsn2
  z = z + dz/2;
  phizx = atan(z/Zrx); phizy = atan(z/Zry);
  scale_z = sqrt(sqrt(1+(z/Zrx)^2)*sqrt(1+(z/Zry)^2));
  Czee = 1i*C0zee.*exp(+1i*(theta_p+dk*z+E0xee*phizx + E0yee*phizy))./scale_z;
  Czeo = 1i*C0zeo.*exp(+1i*(theta_p+dk*z+E0xeo*phizx + E0yeo*phizy))./scale_z;
  Czoe = 1i*C0zoe.*exp(+1i*(theta_p+dk*z+E0xoe*phizx + E0yoe*phizy))./scale_z;
  Czoo = 1i*C0zoo.*exp(+1i*(theta_p+dk*z+E0xoo*phizx + E0yoo*phizy))./scale_z;

  QDsn2ee = Czee*conj((QDsee + QDsn1ee*dz)./2.0);
  QDsn2eo = Czeo*conj((QDseo + QDsn1eo*dz)./2.0);
  QDsn2oe = Czoe*conj((QDsoe + QDsn1oe*dz)./2.0);
  QDsn2oo = Czoo*conj((QDsoo + QDsn1oo*dz)./2.0);
  
  %% calculate Dsn3
  QDsn3ee = Czee*conj((QDsee + QDsn2ee*dz)./2.0);
  QDsn3eo = Czeo*conj((QDseo + QDsn2eo*dz)./2.0);
  QDsn3oe = Czoe*conj((QDsoe + QDsn2oe*dz)./2.0);
  QDsn3oo = Czoo*conj((QDsoo + QDsn2oo*dz)./2.0);
  
  %% calculate Dsn4
  z = z + dz/2;
  phizx = atan(z/Zrx); phizy = atan(z/Zry);
  scale_z = sqrt(sqrt(1+(z/Zrx)^2)*sqrt(1+(z/Zry)^2));

  Czee = 1i*C0zee.*exp(+1i*(theta_p+dk*z+E0xee*phizx + E0yee*phizy))./scale_z;
  Czeo = 1i*C0zeo.*exp(+1i*(theta_p+dk*z+E0xeo*phizx + E0yeo*phizy))./scale_z;
  Czoe = 1i*C0zoe.*exp(+1i*(theta_p+dk*z+E0xoe*phizx + E0yoe*phizy))./scale_z;
  Czoo = 1i*C0zoo.*exp(+1i*(theta_p+dk*z+E0xoo*phizx + E0yoo*phizy))./scale_z;

  QDsn4ee = Czee*conj((QDsee + 2*dz*QDsn3ee))/2;
  QDsn4eo = Czeo*conj((QDseo + 2*dz*QDsn3eo))/2;
  QDsn4oe = Czoe*conj((QDsoe + 2*dz*QDsn3oe))/2;
  QDsn4oo = Czoo*conj((QDsoo + 2*dz*QDsn3oo))/2;

  QDseo = QDseo + dz/3.0.*(QDsn1eo + 2.0*QDsn2eo + 2.0*QDsn3eo + QDsn4eo);
  QDsee = QDsee + dz/3.0.*(QDsn1ee + 2.0*QDsn2ee + 2.0*QDsn3ee + QDsn4ee);
  QDsoe = QDsoe + dz/3.0.*(QDsn1oe + 2.0*QDsn2oe + 2.0*QDsn3oe + QDsn4oe);
  QDsoo = QDsoo + dz/3.0.*(QDsn1oo + 2.0*QDsn2oo + 2.0*QDsn3oo + QDsn4oo);
  
  zpts = zpts + 1
  
end
clear Czee Czoo Czoe Czeo SZ2
clear C0zee C0zeo E0xee E0xoe E0yee E0yeo QDsn1ee QDsn1eo
clear QDsn2ee QDsn2eo QDsn3eo QDsn3ee QDsn4ee QDsn4eo % QDse{e,o} is left
clear C0zoo C0zoe E0xoo E0xeo E0yoo E0yoe QDsn1oo QDsn1oe
clear QDsn2oo QDsn2oe QDsn3oo QDsn3oe QDsn4oo QDsn4oe % QDso{e,o} is left

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

%% calculate the eigenmodes 'v', and squeezing spectrum 'l'.
QDs = [real(QDs); imag(QDs)];
G = full(QDs*(QDs.'));
opts.isreal = 1; opts.issym = 1;
%% [evec,l] = eig(G);
if ( Mx + My < 75 )
	[evec,l] = eig(G);
else
	[evec,l] = eigs( G,100,'sm',opts );
end
l = diag(l);
time_lapse = toc

flipud(l)
if ( dosave )
  save('-mat',sprintf('MV11_Eigenvalue_q%02d_m%02d_Mx%02d_My%02d-P0%g-A0x%g-A0y%g-%s.mat',q,m,Mx,My,P0,ap0x/1e-6,ap0y/1e-6,date),'l','time_lapse','evec');
  save('-mat',sprintf('MV11_GreensFn_q%02d_m%02d_Mx%02d_My%02d-P0%g-A0x%g-A0y%g-%s.mat',q,m,Mx,My,P0,ap0x/1e-6,ap0y/1e-6,date))
end

