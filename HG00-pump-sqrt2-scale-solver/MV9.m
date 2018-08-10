% 
% (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
%     Nonlinear Optics & Nanophotonics Lab, UT-Arlington.
%            
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
%            Pump power is critical value in region of Po = 1e5-1e6.
%            
% 10/29/09 : Remove loop invariants outside the ODE stepper.
%            The transfer matrix evolution is independent of initial
%            conditions which allows us to solve all 2*Mx*My initial conditions
%            at once for obtaining the Greens functions. Fast solver gives
%            Greens function itself in a time of 14.65s for Mx,y=7!!.
%            
%
% 11/04/09 : Modify PSA parameters for NU-parameters:
%            1560nm signal, 780nm pump, 1kw Pump, power, and 20mm
%            crystal with 100x100micron^2 spot size. Prof. V says earlier gain
%            proportional to pump power, and crystal length was quite high
%            which makes the calculations accumulate error, and give 
%            unreasonable squeezing. Also gain must be capped physically
%            realizable to about 30dB or about 1-1000, factor.
%            Reducing pump power, eliminates using RK-45 and Predictor-Corrector
%            methods for the solutions, as relative convergence of step size
%            is eliminated. Eigenvalue spectrum has 1/lambda, lambda dependence
%            between squeezed and anti-squeezed modes with excellent agreement.
% 
% 11/09/09 : Add correct form of phase mismatch factor in evolution.
%            which is independent of 'dz', and just a function of 'z'.
% 
% 11/13/09 : Completely redo the indexing notations. Fix more index
%            notations while calculating DPR, and DLN.
%            This brings me closer to Mathematica output, and Mathematica
%            output closer to this one. Convergence.
% 
% 11/18/09 : Convert to Matlab compatible scripts.
% 
% 11/29/09 : Script corrections from Prof. Vasilyev.
% 12/01/09 : Check differences with Prof. V scripts.
% 12/02/09 : Use sparse matrices at most places. Remove 50% memory hit.
%            Copy some of Prof. Vasilyev design choices to keep C0 real,
%            and move complex parameters inside ODE-IVP loop. Remove
%            conj_op. Generally optimize for maximum memory lean-ness.
%            Clear variables beyond their usage points.
%            But dilemma is to see if adding sparse-ness i/home/muthuncreases
%            memory reach or does it steepen the calculation times?
%            Replace extraneous variables with minimum set.
%	     Remove superfluous QDsn = 0*QDsn type variable reset in loops.
% 	     Remove PLRN
% 12/07/09 : Use complex number matrix ops which are memory efficient.
%            Factor out the odd & even matrices on Cz for speed in matrix
%            products Cz*QDs as there are 50% zeros in this.
%            Use chequered-block matrices.
%            16x16 mode solves in about 104s = 1min 44s! Compare to previous solver
%            at 196s = 3min 16s. And then factor in the 1/2 memory
%            reduction./home/muthu
%
%            Inline the blockset/get scripts.
%            Actually another level of abstraction is possible,
%            by taking alternative rows & columns of the QDse
%            and the QDso. Two level decomposition first by My x My blocks
%            then by odd-even rows gives the optimal performance of the solver.
%            We solve 16x16 mode case in 51s!
% 
% 12/08/09 : Use kronecker products toeplitz and hankel matrix 
%            forms to generate E0x{o,e}, E0y{o,e} respectively.
%            Reorganize to use least memory while downsampling 
%            to calculate C0z{e,o} before the RK-iteration.
%
% 12/10/09 : Factorial2 function fails for arguments above 85 the way
%            it is implemented. So directly compute the Bmm'(z=0)
%            coupling coefficients using gammaln in log_base_e and
%            invert it back. Not memory efficient, but correct nevertheless.
%            
% 12/11/09 : Worry about Matlab compatibility for HPC.
% 
% 12/14/09 : Sparse matrix donot make sense in QDs{o,e}x{o,e} variables./home/muthu
%            Also speye -> eye, as matrices eventual become full. Only 
%            DLN, DPR variable sparse-ness IS useful even for 50% fill-in.
% 
% 01/14/10 : Study a 512x32 case.
% 
% 02/18/10 : U9.9e-12;se eigs and just study the partial eigenspectra starting from the
%            most squeezed mode to the successive 50-60 eigenvalues.
% 
% 02/18/10 : Choose the power where the gain = 15, for 400x100 micron^2
%            spot size. Generate cases to check with the older data
%            from Prof. V, and validate MV9, MV9_parallel solvers.
% 03/07/10 :  Make G from sparse -> full matrix.
%
% 03/24/10 : The problem is with using sparse matrix input to blockset(); this
%            is terribly slow. We need to use full matrices input at this
%            point. Also eigs() works just fine for limited number of
%            eigenvectors. Make ODE loop matrices { C0z{e}x{o} }
%            full-matrix.
% 
% 06/05/10 : Add the 'dosave' variable for specific data saving runs.
%            Save full eigenvalue spectrum for small cases whenever youre calculating them.
%            All just like MV10.m.
% 
% 11/10/10 : Use xi_eff = 2*deff; and use Kappa with xi_eff.
% 11/16/10 : Add deff to the filenames for consistency.
% 
% 12/29/10 : Run 100x100 spot size case at dk = 220m^-1
%            at 5,7,9,10,12,15,16 kW powers.
% 
clear all
close all
more off
echo off

% start after a few hours to ease the load.
% pause(12*3600);

tic
%% Number of signal-idler modes.

%% Mx = 16; My = 16;
%Mx = 128; My = 128;
%%Mx = 4; My = 4;
%%Mx = 32; My = 32;
%Mx = 16; My = 16;
%Mx = 6; My = 4;
%Mx = 44; My = 4;
%Mx = 2; My = 2;
Mx = 32; My = 32;
%%Mx = 64; My = 32;
%%Mx = 1024; My = 32;
%Mx = 512; My = 32;
%Mx = 128; My = 32;
%Mx = 4; My = 2;
%Mx = 128; My = 128; %for 200x200 spot
%%Mx = 2048; My = 8; %for 800x50 spot
%%Mx = 128; My = 32;

% for 400x25, 400x50 or 350x25, 350x50
%Mx = 512; My = 32;

% for 300x25 spot 
%Mx = 384; My = 32;

dosave = 1;
if ( mod(Mx,2) + mod(My,2) ~= 0 )
  error(' Mx, My need to be even simultaneously ')
end

%% using MKS units
c=2.99792458e8; eps0=8.85e-12;

%% NU-experiment parameters
%% Lz : 20mm crystal, signal: 1560nm, pump: 780nm, deff : 8.7e-12, ns,p: 1.78
lambda_s = 1560e-9; lambda_p = lambda_s/2; ks = 2*pi/lambda_s;  kp = ks*2;
%data for KTP
%Lz =  20e-3; ns = 1.78; np = ns; deff = 8.7e-12;
%xi_eff = 2*deff; omega_s = ks*c;

% Type-0, 3cm PPKTP Gideon
fprintf('Type-0 3cm PPKTP / Gideon\n')
Lz = 30e-3; ns = 1.78; np = ns;
% [5.7,6,6.3]
deff = 5.7e-12; xi_eff = 2*deff; omega_s = ks*c;
P0 = 0.25e3; ap0x=27e-6; ap0y=25.7e-6;

% Type-II PPKTP / Gideon
%fprintf('PPKTP / Gideon\n')
%Lz =  20e-3; ns = 1.78; np = ns; 
%deff = 2.5e-12; xi_eff = 2*deff; omega_s = ks*c;

%data for PPLN (pip-lin) / Amar
%fprintf('PPLN / Amar\n')
%Lz = 20e-3; ns = 2.14; np = ns; 2.19; 
%deff = 8.7e-12; xi_eff = 2*deff; omega_s = ks*c;

%% BBN parameters for the 400x25, 200x25, 400x50, 200x50
%% spot sizes, and also current 800x25.
%Lz =  20e-3; ns = 1.78; np = ns; deff = 8.7e-12; xi_eff = 2*deff; omega_s = ks*c;

%% Harris Experiments on 07/18/11
%Lz = 25e-3; ns = 2.14; np = ns; deff = 8.0e-12; xi_eff = 2*deff; omega_s = ks*c;

dk = 4*pi*(np-ns)/lambda_s;
fprintf('Presetting dk value ...\n');
%dk = -110; -220;+220; %200m^-1
fprintf('Setting dk = %g\n',dk);

theta_p = -pi/2;

%% simulation parameters
Lpts = 300; dz = Lz/(Lpts);

%P0 = 10e3; %5.4e3; 
%ap0x = 800e-6; ap0y = 25e-6;

%% Harris Experiments on 07/18/11
%P0 = 0.185e3; ap0x=19.5e-6; ap0y=ap0x;

%% pump parameters: NU PPKTP setup
%P0 = 0.25*1e3; 1e3;  ap0x=36.8e-6; ap0y=36.8e-6; %ap0x = 40e-6; ap0y = ap0x;

%% pump parameters: NU PPLN setup;
%P0 = 0.625e3; 3e3; ap0x = 39e-6; ap0y = 39e-6;

%P0 = 4e3; ap0x = 200e-6; ap0y=200e-6;
%P0 = 1.25e3;ap0x = 150e-6; ap0y = 100e-6;
%for P0 = linspace(5.25,10.5,5)*1e3 % kW powers.

%% pump parameters: 20kW pump case with 1/e waist size of pump
%%P0 = 17.25e3; ap0x = 400e-6; ap0y = 100e-6;
%%P0 = 16.5e3; ap0x = 200e-6; ap0y = 200e-6;
%%P0 = 10.00e3; ap0x = 200e-6; ap0y = 100e-6;
%%P0 = 16.00e3; ap0x = 800e-6; ap0y = 50e-6;
%%P0 = 8.75e3; ap0x = 200e-6; ap0y = 100e-6;
%%P0 = 9.00e3; ap0x = 200e-6; ap0y = 100e-6;
%%P0 = 19.50e3; ap0x = 800e-6; ap0y = 50e-6;
%%P0 = 16e3; ap0x = 800e-6; ap0y = 50e-6;

%% Data points for the Harris setup
%P0 = 3e3; ap0x = 200e-6; ap0y = 25e-6;
%P0 = 3e3; ap0x = 350e-6; ap0y = 25e-6;

%% 08/11/10 compare HO data @ smaller pump with the fundamental pump
%%P0 = 16.25e3; ap0x = 700e-6; ap0y = 50e-6;

%P0 = 2.4e3; ap0x = 100e-6; ap0y = 100e-6;

%% Calibration with data for CLEO-10, and PQE-10.
%% Powers adjusted for uniform gain of 15.
%%
%P0 = (1.15e3)/4; ap0x = 25e-6; ap0y = 25e-6;
%%P0 = 2.4e3; ap0x = 100e-6; ap0y = 25e-6;
%%P0 = 3.1e3; ap0x = 100e-6; ap0y = 50e-6;
%%P0 = 5e3; ap0x = 100e-6; ap0y = 100e-6;
%%P0 = 16.25e3; ap0x = 200e-6; ap0y = 200e-6;
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
fprintf('Running case : P0 = %g [kW], spot : %d x %d at modes %dx%d \n',P0/1e3, ap0x/1e-6,ap0y/1e-6,Mx,My);
fprintf('d_eff = %g pm/V \n',deff/1e-12);
fprintf('Fundamental Order Pump mode\n');
fprintf('*************************************************************\n')

as0x = sqrt(2)*ap0x; as0y = sqrt(2)*ap0y;

%% Rayleigh length for X, Y modes
Zrx = 2*pi*np*ap0x^2/lambda_p; Zry = 2*pi*np*ap0y^2/lambda_p;
Kappa = omega_s*xi_eff/sqrt(2*eps0*ns^2*np*c^3)*sqrt(P0/(pi*ap0x*ap0y));

%% fast Gouy phase calculation 
evenMx = 0.5 + (0:2:2*Mx);
E0xe = kron(hankel(evenMx(1:Mx/2),evenMx(Mx/2:Mx-1)),toeplitz(mod(1:My,2)));
E0xo = kron(2 + hankel(evenMx(1:Mx/2),evenMx(Mx/2:Mx-1)),toeplitz(mod(1:My,2)));
oddMy = (0.5 + (0:(2*My-1))).*mod([1:My*2],2);
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
LN = kronsum( lL', nN );
nzLN = sparse(mod(LN,2)==0);
LL = ( LL.*nzLN ); LN = ( LN.*nzLN );
DLN = sparse(real((-1).^((LL-NN)/2)).*exp(gammaln(LN + 1) - ...
    (  gammaln(1+LN/2) + LN*log(2) + 0.5*(log(2)+gammaln( 1 + LL ) + ...
					  gammaln( 1 + NN )))));
DLN(kron(toeplitz(mod(0:Mx-1,2)),ones(My))>0)=0;
[DLNe,DLNo] = blockget(  DLN, Mx, My );
DLNee = DLNe(1:2:end,1:2:end); DLNeo = DLNe(2:2:end,2:2:end);
DLNoe = DLNo(1:2:end,1:2:end); DLNoo = DLNo(2:2:end,2:2:end);
clear NN LL LN nzLN DLN DLNe DLNo

PP = repmat( pP', [1, Mx*My] ); RR = PP';
PR = kronsum( pP', rR );
nzPR = sparse( mod(PR,2)==0);
PP = ( PP.*nzPR ); PR = ( PR.*nzPR );
DPR = sparse(real((-1).^((PP-RR)/2)).*exp(gammaln(PR + 1) - ...
    (  gammaln(1+PR/2) + PR*log(2) + 0.5*(log(2)+gammaln( 1 + PP ) + gammaln( 1 + RR )))));
DPR(kron(ones(Mx),toeplitz(mod([0:(My-1)],2)))>0)=0;
[DPRe,DPRo] = blockget(  DPR, Mx, My );
DPRee = DPRe(1:2:end,1:2:end); DPReo = DPRe(2:2:end,2:2:end);
DPRoe = DPRo(1:2:end,1:2:end); DPRoo = DPRo(2:2:end,2:2:end);
clear RR PP PR nzPR DPR DPRe DPRo
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
zz = [z];
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
  zz = [ zz; z];
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
  zz=[zz;z];
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
format long
fprintf('%g\n',zz)
%echo on
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
save('-mat',sprintf('MV9_Eigenvalue_Mx%02d_My%02d-deff-%g-P0%g-A0x%g-A0y%g-%s.mat',Mx,My,deff/1e-12,P0,ap0x/1e-6,ap0y/1e-6,date),'l','time_lapse','evec');
save('-v7.3',sprintf('MV9_GreensFn_Mx%02d_My%02d-deff-%g-P0%g-A0x%g-A0y%g-%s.mat',Mx,My,deff/1e-12,P0,ap0x/1e-6,ap0y/1e-6,date))
end
%end
