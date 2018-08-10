% 
% (C) 2009 Muthiah Annamalai <muthiah.annamalai@mavs.uta.edu>
%     Nonlinear Optics & Nanophotonics Lab, UT-Arlington.
%            
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% 10/27/09 : PSA modal theory solver for DOPA as formulated in memo 
%            by Prof. Vasilyev. Prof. V's new version of overlap integrals
%            includes the selection rules automatically.
%            
%            DOPA LG mode solver based on closed-form expression
%            for a LG00 pump. The Kumar-Vasilyev conventions are
%            used for the field dependence on +z, and frequency,
%            while the waist is located at 1/e intensity on transverse
%            plane.
%            
%            Convert to using quadratures of the modes.
%            
%            Calculate the Green's functions.
%            
% 05/11/10 : We modify MV9 solver to use as0x = k*ap0x, when k != sqrt(2)
%            as previously assumed, requiring recalculation of the overlap
%            integrals at each space point.
% 
% 05/19/10 : Adding bones of solver. Breathe in the life, finally.
%            Works but definitely wrong; need to check the reorganizations.
% 
% 06/3/10  : Fix some programming errors with correct formulas for overlap integral.
% 06/4/10  : Add new variable dosave. We try to get reasonable number of
%            eigenvalues & eigenmodes or full spectrum of eigenvalues
%            depending on problem size.
% 06/09/11 : I try to fix the overflow problem, by calculating the
%            overlap integrals correctly like in
%            overlap_sqrt2_unscaled.m script, using new scaled
%            coefficients imported from Mathematica.
% 
%            Move to newer values of pump powers, which is xeff = 2*deff
% 
% 09/13/11 : I have identified some small but imporant corrections to the
%            overlap integral formulas. This may just make this script work
%            finally!
% 
clear all
close all
more off
echo off
tic
%% Number of signal-idler modes.

%% Mx = 16; My = 16;
%% Mx = 4; My = 4;
%Mx = 32; My = 32;
%Mx = 16; My = 16;
%Mx = 6; My = 4;
%Mx = 4; My = 6;
%Mx = 2; My = 4;
Mx = 32; My = 32;
%%Mx = 64; My = 32;
%%Mx = 1024; My = 32;
%Mx = 512; My = 32;
% Mx = 128; My = 128;
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
Lz =  20e-3; ns = 1.78; np = ns; deff = 8.7e-12; xeff = 2*deff; omega_s = ks*c;
dk = 4*pi*(np-ns)/lambda_s;
theta_p = -pi/2;

%% simulation parameters
Lpts = 300; dz = Lz/(Lpts);

%% Calibration with data for CLEO-10, and PQE-10.
%% Powers adjusted for uniform gain of 15.
P0 = (1.15e3)/4; ap0x = 25e-6; ap0y = 25e-6;
%%P0 = 2.4e3; ap0x = 100e-6; ap0y = 25e-6;
%%P0 = 3.1e3; ap0x = 100e-6; ap0y = 50e-6;
%%P0 = 5e3; ap0x = 100e-6; ap0y = 100e-6;
%%P0 = 16.25e3; ap0x = 200e-6; ap0y = 200e-6;
%%P0 = 21.5e3; ap0x = 800e-6; ap0y = 50e-6;
%%P0 = 17.25e3; ap0x = 400e-6; ap0y = 100e-6;

%kscalex = sqrt(2); kscaley = sqrt(2);
%as0x = kscalex*ap0x; as0y = kscaley*ap0y;

%% optimized signal waist
as0x = ap0x*sqrt(2);%can compare with MV9 at f_s = Sqrt[2] case
as0x = 20e-6; 32.1572e-6; 31.4e-6; as0y = as0x; %compact mode calculated from MV9
kscalex = as0x/ap0x; kscaley = as0y/ap0y;
fprintf('scale factors are (x,y) = %g \n',[kscalex,kscaley]);

%% Rayleigh length for X, Y modes
Zrpx0 = 2*pi*np*ap0x^2/lambda_p; Zrsx0 = (kscalex)^2*Zrpx0/2; 
Zrpy0 = 2*pi*np*ap0y^2/lambda_p; Zrsy0 = (kscaley)^2*Zrpy0/2;
Kappa = omega_s*xeff/sqrt(2*eps0*ns^2*np*c^3)*sqrt(P0);

%% 
%% OVERLAP INTEGRAL CALCULATION
%%
%% CURRENT METHOD
%% Data values of scaled Hermite polynomial convolution
%% coefficients from order 0-32 with itself.
%% Pre-scaled by value 1/(pi^(3/4)*sqrt(2^(p+r)*factorial(p)*factorial(r)))
%% Adjust the scale by 1/sqrt(ap0x) for every specific spot size.
%% OLDER METHOD
%% overlap integral calculation Lookup Table including the scale factor
%% of 1/(pi*sqrt(pi)*2^(p+r)*r!*p!) on each discrete convolved coefficient.
%MMX = 64; MMY = 64;
%HermiCoeffs = load(sprintf('ConvCoeffHermiScaled_%d_%d.dat',MMX,MMY));

MMX=33; MMY=33;
load('UnscaledOverlap_Coeffs_Mx33.mat');
HermiCoeffs=UnscaledOverlap_Coeffs_Mx33;
scale_x = 1/sqrt(ap0x); scale_y = 1/sqrt(ap0y);

fprintf('Entering re-indexing loop\n')
for i = 0:Mx*My-1;
    p = mod( i, My ); l = floor( i / My );
    for j = 0:Mx*My-1;
        r = mod( j, My ); n = floor( j / My );
        polyidxPR(1+i,1+j) = p*MMY + r + 1;
        polyidxLN(1+i,1+j) = l*MMY + n +1;
        EPR(1+i,1+j) = p + r + 1;
	ELN(1+i,1+j) = l + n + 1;
    end
end
clear i j p r l n
W = -0.5-[0:2:(Mx+My)]/2;
calc_overlap_x = @( xix ) scale_x*[(xix.^(W))]*[(HermiCoeffs(polyidxPR,1:2:(Mx+My+1))).'];
calc_overlap_y = @( xiy ) scale_y*[(xiy.^(W))]*[(HermiCoeffs(polyidxLN,1:2:(Mx+My+1))).'];

%% fast Gouy phase calculation 
evenMx = 1 + (0:2:2*Mx);
E0xe = kron(hankel(evenMx(1:Mx/2),evenMx(Mx/2:Mx-1)),toeplitz(mod(1:My,2)));
E0xo = kron(2 + hankel(evenMx(1:Mx/2),evenMx(Mx/2:Mx-1)),toeplitz(mod(1:My,2)));
oddMy = (1 + (0:(2*My-1))).*mod([1:My*2],2);
E0ye = kron( ones(Mx/2,Mx/2), hankel( oddMy(1:My), oddMy(My:2*My-1)));
E0yo = E0ye;
 
E0xee = E0xe(1:2:end,1:2:end); E0xeo = E0xe(2:2:end,2:2:end);
E0xoe = E0xo(1:2:end,1:2:end); E0xoo = E0xo(2:2:end,2:2:end);

E0yee = E0ye(1:2:end,1:2:end); E0yeo = E0ye(2:2:end,2:2:end);
E0yoe = E0yo(1:2:end,1:2:end); E0yoo = E0yo(2:2:end,2:2:end);

clear E0ye E0yo E0xo E0xe evenMx oddMy

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

% Indexing Conventions for 4D ((p,l),(r,n)) -> 2D (i,j) contraction
% notation. In new index conventions, i is row index, and j is col index,
% which represent how the j-th mode couples to i-th mode.
%
% Row : i == (p,l) mode  p : row index, l : col index
%       p = mod( i, My ); l = floor( i / My );
% Col : j == (r,n) modsparsee  transpose of ( r : row index, n : col index )
%       r = mod( j, My ); n = floor( j / My );
% keep everything row-vectors
%
z = -Lz/2; zpts = 0;
fprintf('Entering z-iteration loop\n')
while ( z < Lz/2 & zpts < Lpts )
    fprintf('z step [ %d] at distance %g mm\n',zpts,z/1e-3);

    [max(max(abs(QDseo))), max(max(abs(QDsee))), max(max(abs(QDsoe))), max(max(abs(QDsoo)))]
    
    %% Coupling integrals become Z-dependent.    
    aszx = as0x*sqrt(1 + (z/Zrsx0)^2); apzx = ap0x*sqrt(1 + (z/Zrpx0)^2);
    if ( abs(z) < 10*eps )
        xi_x = (1+(1/2 )*(aszx/apzx)^2 );
    else
        Rpzx = z*(1+(Zrpx0/z)^2);
        Rszx = z*(1+(Zrsx0/z)^2);
        xi_x = 1+((1/2)*(aszx/apzx)^2 ) ...
            - 1i*( ks*as0x^2*(1+(z/Zrsx0)^2)*(1/Rpzx - 1/Rszx) );
    end
    
    % calculate all the overlap coefficients
    %%DPR = [(xi_x.^(W)).*gamma(-W)]*[(HermiCoeffs(polyidxPR,1:2:(Mx+My+1)))'];
    %DPR = [(xi_x.^(W))]*[(HermiCoeffs(polyidxPR,1:2:(Mx+My+1))).'];
    DPR = calc_overlap_x( xi_x);
    max(abs(DPR))
    DPR = reshape( DPR, Mx*My,Mx*My);
    max(max(abs(DPR)))
    %DPR(kron(ones(Mx),toeplitz(mod([0:(My-1)],2)))>0)=0;
    [DPRe,DPRo] = blockget(  DPR, Mx, My );
    DPRee = DPRe(1:2:end,1:2:end); DPReo = DPRe(2:2:end,2:2:end);
    DPRoe = DPRo(1:2:end,1:2:end); DPRoo = DPRo(2:2:end,2:2:end);
    clear DPR DPRe DPRo
    
    aszy = as0y*sqrt(1 + (z/Zrsy0)^2); apzy = ap0y*sqrt(1 + (z/Zrpy0)^2);
    if ( abs(z) < 10*eps )
        xi_y = (1+(1/2)*(aszy/apzy)^2 );
    else
        Rpzy = z*(1+(Zrpy0/z)^2);
        Rszy = z*(1+(Zrsy0/z)^2);
        xi_y = (1+(1/2)*(aszy/apzy)^2 ) ...
            - 1i*( ks*as0y^2*(1+(z/Zrsy0)^2)*(1/Rpzy - 1/Rszy) );
    end
    
    % calculate all the overlap coefficients
    %DLN = [(xi_y.^(W)).*gamma(-W)]*[(HermiCoeffs(polyidxLN,1:2:(Mx+My+1)))'];
    DLN = calc_overlap_y(xi_y);
    max(abs(DLN))
    DLN = reshape( DLN, Mx*My,Mx*My);
 
    % DLN(kron(toeplitz(mod(0:Mx-1,2)),ones(My))>0)=0;
    [DLNe,DLNo] = blockget(  DLN, Mx, My );
    DLNee = DLNe(1:2:end,1:2:end); DLNeo = DLNe(2:2:end,2:2:end);
    DLNoe = DLNo(1:2:end,1:2:end); DLNoo = DLNo(2:2:end,2:2:end);
    clear DLN DLNe DLNo
    
    C0zee = full(Kappa*DPRee.*DLNee);
    C0zoe = full(Kappa*DPRoe.*DLNoe);
    C0zeo = full(Kappa*DPReo.*DLNeo);
    C0zoo = full(Kappa*DPRoo.*DLNoo);

%    [max(max(abs(C0zeo))), max(max(abs(C0zee))), max(max(abs(C0zoe))), max(max(abs(C0zoo)))]
    clear DPRoe DPRoo DPRee DPReo
    clear DLNoo DLNoe DLNeo DLNee
    %end of overlap integral calculation for this z-step.
    
    %% calculate Dsn1
    phizsx = atan(z/Zrsx0); phizsy = atan(z/Zrsy0);
    phizpx = atan(z/Zrpx0); phizpy = atan(z/Zrpy0);
    scale_z = sqrt(sqrt(1+(z/Zrpx0)^2)*sqrt(1+(z/Zrpy0)^2));
    Czee = 1i*C0zee.*exp(+1i*(theta_p+dk*z+E0xee*phizsx - phizpx/2 - phizpy/2 + E0yee*phizsy))./scale_z;
    Czeo = 1i*C0zeo.*exp(+1i*(theta_p+dk*z+E0xeo*phizsx  - phizpx/2 - phizpy/2 + E0yeo*phizsy))./scale_z;
    Czoe = 1i*C0zoe.*exp(+1i*(theta_p+dk*z+E0xoe*phizsx  - phizpx/2 - phizpy/2 + E0yoe*phizsy))./scale_z;
    Czoo = 1i*C0zoo.*exp(+1i*(theta_p+dk*z+E0xoo*phizsx  - phizpx/2 - phizpy/2 + E0yoo*phizsy))./scale_z;
    
    QDsn1ee = Czee*conj(QDsee)./2.0;
    QDsn1eo = Czeo*conj(QDseo)./2.0;
    QDsn1oe = Czoe*conj(QDsoe)./2.0;
    QDsn1oo = Czoo*conj(QDsoo)./2.0;
    
    [max(max(abs(QDseo))), max(max(abs(QDsee))), max(max(abs(QDsoe))), max(max(abs(QDsoo)))]        
    pause(1)	     

    %% calculate Dsn2
    z = z + dz/2;
    
    %% Coupling integrals become Z-dependent.    
    aszx = as0x*sqrt(1 + (z/Zrsx0)^2); apzx = ap0x*sqrt(1 + (z/Zrpx0)^2);
    if ( abs(z) < 10*eps )
        xi_x = (1+(1/2)*(aszx/apzx)^2 );
    else
        Rpzx = z*(1+(Zrpx0/z)^2);
        Rszx = z*(1+(Zrsx0/z)^2);
        xi_x = (1+(1/2)*(aszx/apzx)^2 ) ...
            - 1i*( ks*as0x^2*(1+(z/Zrsx0)^2)*(1/Rpzx - 1/Rszx) );
    end
    
    % calculate all the overlap coefficients
    %%DPR = [(xi_x.^(W)).*gamma(-W)]*[(HermiCoeffs(polyidxPR,1:2:(Mx+My+1)))'];
    %%DPR = [(xi_x.^(W))]*[(HermiCoeffs(polyidxPR,1:2:(Mx+My+1))).'];
    DPR = calc_overlap_x( xi_x);
    max(abs(DPR))
    DPR = reshape( DPR, Mx*My,Mx*My);
    
    %DPR(kron(ones(Mx),toeplitz(mod([0:(My-1)],2)))>0)=0;
    [DPRe,DPRo] = blockget(  DPR, Mx, My );
    DPRee = DPRe(1:2:end,1:2:end); DPReo = DPRe(2:2:end,2:2:end);
    DPRoe = DPRo(1:2:end,1:2:end); DPRoo = DPRo(2:2:end,2:2:end);
    clear  DPR DPRe DPRo  
    
    aszy = as0y*sqrt(1 + (z/Zrsy0)^2); apzy = ap0y*sqrt(1 + (z/Zrpy0)^2);
    if ( abs(z) < 10*eps )
        xi_y = (1+(1/2)*(aszy/apzy)^2 );
    else
        Rpzy = z*(1+(Zrpy0/z)^2);
        Rszy = z*(1+(Zrsy0/z)^2);
        xi_y = (1+(1/2)*(aszy/apzy)^2 ) ...
            - 1i*( ks*as0y^2*(1+(z/Zrsy0)^2)*(1/Rpzy - 1/Rszy) );
    end
    
    % calculate all the overlap coefficients
    %DLN = [(xi_y.^(W)).*gamma(-W)]*[(HermiCoeffs(polyidxLN,1:2:(Mx+My+1)))'];
    DLN = calc_overlap_y(xi_y);
    max(abs(DLN))
    DLN = reshape( DLN, Mx*My,Mx*My);
    
    %DLN(kron(toeplitz(mod(0:Mx-1,2)),ones(My))>0)=0;
    [DLNe,DLNo] = blockget(  DLN, Mx, My );
    DLNee = DLNe(1:2:end,1:2:end); DLNeo = DLNe(2:2:end,2:2:end);
    DLNoe = DLNo(1:2:end,1:2:end); DLNoo = DLNo(2:2:end,2:2:end);
    clear DLN DLNe DLNo
    
    C0zee = full(Kappa*DPRee.*DLNee);
    C0zoe = full(Kappa*DPRoe.*DLNoe);
    C0zeo = full(Kappa*DPReo.*DLNeo);
    C0zoo = full(Kappa*DPRoo.*DLNoo);
    
    clear DPRoe DPRoo DPRee DPReo
    clear DLNoo DLNoe DLNeo DLNee
    %end of overlap integral calculation for this z-step.
    
    phizsx = atan(z/Zrsx0); phizsy = atan(z/Zrsy0);
    phizpx = atan(z/Zrpx0); phizpy = atan(z/Zrpy0);
    scale_z = sqrt(sqrt(1+(z/Zrpx0)^2)*sqrt(1+(z/Zrpy0)^2));
    
    Czee = 1i*C0zee.*exp(+1i*(theta_p+dk*z+E0xee*phizsx  - phizpx/2 - phizpy/2 + E0yee*phizsy))./scale_z;
    Czeo = 1i*C0zeo.*exp(+1i*(theta_p+dk*z+E0xeo*phizsx  - phizpx/2 - phizpy/2 + E0yeo*phizsy))./scale_z;
    Czoe = 1i*C0zoe.*exp(+1i*(theta_p+dk*z+E0xoe*phizsx  - phizpx/2 - phizpy/2 + E0yoe*phizsy))./scale_z;
    Czoo = 1i*C0zoo.*exp(+1i*(theta_p+dk*z+E0xoo*phizsx  - phizpx/2 - phizpy/2 + E0yoo*phizsy))./scale_z;
    
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
    
    %% Coupling integrals become Z-dependent.    
    aszx = as0x*sqrt(1 + (z/Zrsx0)^2); apzx = ap0x*sqrt(1 + (z/Zrpx0)^2);
    if ( abs(z) < 10*eps )
        xi_x = (1+(1/2 )*(aszx/apzx)^2 );
    else
        Rpzx = z*(1+(Zrpx0/z)^2);
        Rszx = z*(1+(Zrsx0/z)^2);
        xi_x = (1+(1/2)*(aszx/apzx)^2 ) ...
            - 1i*( ks*as0x^2*(1+(z/Zrsx0)^2)*(1/Rpzx - 1/Rszx) );
    end
    
    % calculate all the overlap coefficients
    %DPR = [(xi_x.^(W)).*gamma(-W)]*[(HermiCoeffs(polyidxPR,1:2:(Mx+My+1)))'];
    %%DPR = [(xi_x.^(W))]*[(HermiCoeffs(polyidxPR,1:2:(Mx+My+1))).'];
    DPR = calc_overlap_x( xi_x);
    max(abs(DPR))
    DPR = reshape( DPR, Mx*My,Mx*My);
    
    %DPR(kron(ones(Mx),toeplitz(mod([0:(My-1)],2)))>0)=0;
    [DPRe,DPRo] = blockget(  DPR, Mx, My );
    DPRee = DPRe(1:2:end,1:2:end); DPReo = DPRe(2:2:end,2:2:end);
    DPRoe = DPRo(1:2:end,1:2:end); DPRoo = DPRo(2:2:end,2:2:end);
    clear DPR DPRe DPRo
    
    aszy = as0y*sqrt(1 + (z/Zrsy0)^2); apzy = ap0y*sqrt(1 + (z/Zrpy0)^2);
    if ( abs(z) < 10*eps )
        xi_y = (1+(1/2)*(aszy/apzy)^2 );
    else
        Rpzy = z*(1+(Zrpy0/z)^2);
        Rszy = z*(1+(Zrsy0/z)^2);
        xi_y = (1+(1/2)*(aszy/apzy)^2 ) ...
            - 1i*( ks*as0y^2*(1+(z/Zrsy0)^2)*(1/Rpzy - 1/Rszy) );
    end
    
    % calculate all the overlap coefficients
    %DLN = [(xi_y.^(W)).*gamma(-W)]*[(HermiCoeffs(polyidxLN,1:2:(Mx+My+1)))'];
    DLN = calc_overlap_y(xi_y);
    max(abs(DLN))
    DLN = reshape( DLN, Mx*My,Mx*My);
    
    %DLN(kron(toeplitz(mod(0:Mx-1,2)),ones(My))>0)=0;
    [DLNe,DLNo] = blockget(  DLN, Mx, My );
    DLNee = DLNe(1:2:end,1:2:end); DLNeo = DLNe(2:2:end,2:2:end);
    DLNoe = DLNo(1:2:end,1:2:end); DLNoo = DLNo(2:2:end,2:2:end);
    clear DLN DLNe DLNo
    
    C0zee = full(Kappa*DPRee.*DLNee);
    C0zoe = full(Kappa*DPRoe.*DLNoe);
    C0zeo = full(Kappa*DPReo.*DLNeo);
    C0zoo = full(Kappa*DPRoo.*DLNoo);
    
    clear DPRoe DPRoo DPRee DPReo
    clear DLNoo DLNoe DLNeo DLNee
    %end of overlap integral calculation for this z-step.
    
    phizsx = atan(z/Zrsx0); phizsy = atan(z/Zrsy0);
    phizpx = atan(z/Zrpx0); phizpy = atan(z/Zrpy0);
    scale_z = sqrt(sqrt(1+(z/Zrpx0)^2)*sqrt(1+(z/Zrpy0)^2));
    
    Czee = 1i*C0zee.*exp(+1i*(theta_p+dk*z+E0xee*phizsx  - phizpx/2 - phizpy/2 + E0yee*phizsy))./scale_z;
    Czeo = 1i*C0zeo.*exp(+1i*(theta_p+dk*z+E0xeo*phizsx  - phizpx/2 - phizpy/2 + E0yeo*phizsy))./scale_z;
    Czoe = 1i*C0zoe.*exp(+1i*(theta_p+dk*z+E0xoe*phizsx  - phizpx/2 - phizpy/2 + E0yoe*phizsy))./scale_z;
    Czoo = 1i*C0zoo.*exp(+1i*(theta_p+dk*z+E0xoo*phizsx  - phizpx/2 - phizpy/2 + E0yoo*phizsy))./scale_z;
    
    QDsn4ee = Czee*conj((QDsee + 2*dz*QDsn3ee))/2;
    QDsn4eo = Czeo*conj((QDseo + 2*dz*QDsn3eo))/2;
    QDsn4oe = Czoe*conj((QDsoe + 2*dz*QDsn3oe))/2;
    QDsn4oo = Czoo*conj((QDsoo + 2*dz*QDsn3oo))/2;
    
    [max(max(abs(QDseo))), max(max(abs(QDsee))), max(max(abs(QDsoe))), max(max(abs(QDsoo)))]
    QDseo = QDseo + dz/3.0.*(QDsn1eo + 2.0*QDsn2eo + 2.0*QDsn3eo + QDsn4eo);
    QDsee = QDsee + dz/3.0.*(QDsn1ee + 2.0*QDsn2ee + 2.0*QDsn3ee + QDsn4ee);
    QDsoe = QDsoe + dz/3.0.*(QDsn1oe + 2.0*QDsn2oe + 2.0*QDsn3oe + QDsn4oe);
    QDsoo = QDsoo + dz/3.0.*(QDsn1oo + 2.0*QDsn2oo + 2.0*QDsn3oo + QDsn4oo);
    
    zpts = zpts + 1
    [max(max(abs(QDseo))), max(max(abs(QDsee))), max(max(abs(QDsoe))), max(max(abs(QDsoo)))]
    pause(1)
end
echo on
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
if ( Mx + My < 75 )
   [evec,l] = eig(G);
else
   [evec,l] = eigs( G,100,'sm',opts );
end
%% [evec,l] = eigs( G,10,'sm',opts );
l = diag(l);
time_lapse = toc

flipud(l)
if ( dosave )
       save('-mat',sprintf('MV10_Eigenvalue_Mx%02d_My%02d-kscalex-%g-kscaley-%g-P0%g-A0x%g-A0y%g-%s.mat',Mx,My,kscalex,kscaley,P0,ap0x/1e-6,ap0y/1e-6,date),'l','time_lapse','evec');
       save('-v7.3',sprintf('MV10_GreensFn_Mx%02d_My%02d-kscalex-%g-kscaley-%g-P0%g-A0x%g-A0y%g-%s.mat',Mx,My,kscalex,kscaley,P0,ap0x/1e-6,ap0y/1e-6,date))
end
