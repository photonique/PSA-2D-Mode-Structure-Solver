% (C) 2009 Muthiah Annamalai <muthuspost@gmail.com>
%     Nonlinear Optics & Nanophotonics Lab, UT-Arlington.
%
% 
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
% 09/24/09 : DOPA LG mode solver based on closed-form expression
%            for a LG00 pump. The Kumar-Vasilyev conventions are
%            used for the field dependence on +z, and frequency,
%            while the waist is located at 1/e intensity on transverse
%            plane.
%
% 10/25/09 : d/dz(Dp,l) = sum( [over r,n ] C{p,r,l,n=-l} * Dr,n  )
%            Crn has to be formed using the selection rule n = -l at a
%            given value of p,l.
% 
clear all
close all
% Parameters of calculations
M = 15; % Number of signal-idler modes. Keep it odd.
#M = 3;
M = 5; 
Mn = (M-1);

% all length units in microns
lambda_s = 1064; lambda_p = lambda_s/2; ks = 2*pi/lambda_s;  kp = ks*2;
a0 = 5; %1/e waist size of pump
L =  5.21*1e3; %length of crystal 5.21mm
Lpts = 500; %discretization steps for ODE
Zr = 2*pi*a0^2/lambda_p; %Rayleigh length
dk = 0; ns = 1.5; deff = 8.7e-1;
Ap = sqrt( 1e4 ); %10kW pump
#Ap = 1e-2;

% load the modal coefficients from number3 image.
load number3.mat
[W,H] = size( data );

## s = 1,1.5,2 works out  OK for z = -L/2
## when rmax = ax*s*3, a0 = lambda_s/s;
s = 2; z = -L/2;
a0 = lambda_s/s;  %1/e waist size of pump
Zr = 2*pi*a0^2/lambda_p;
az = a0*sqrt(1 + (z/Zr).^2 );

# rr = linspace(0,az*s*3,W);
# tt = linspace(-pi,pi,H);
# dr = rr(2) - rr(1); dt = tt(2) - tt(1);
# dx = dr*cos( dt ); dy = dr*sin( dt );
# [r,t] = meshgrid( rr, tt );
# X = r.*cos( t );Y = r.*sin( t );

#Ds = flgtx2d( M, M, a0, r, t, z, lambda_s, data,  dx, dy );
Ds = zeros(M,M);
Ds(1,Mn/2+1) = 1;
Ds(1,1) = 10;
fprintf('Initial Signal modal-coeffs\n')
Dsorig = Ds;
abs( Ds )

# orig_data = ilgtx2d( a0, r, t, z, lambda_s, Ds,  dx, dy );
# figure
# imagesc( abs( orig_data ) )
# figure
# imagesc( abs( data ) )
# sum( sum( abs(orig_data - data).^2 ) )
# pause

C0 = zeros(M,M,M);%coupling integral
E0 = zeros(M,M,M);%weight for exponent
f1 = 2/sqrt(pi)*1/a0;
for p = 0:M-1
    for r = 0:M-1
        for l = -Mn/2:+Mn/2
            labs = abs(l);
            f2 = factorial(p+r+labs);
            f3 = sqrt( factorial(p ) * factorial( r) * factorial( p + labs) * factorial( r+ labs));
            f4 = 1/2^(p+r+labs+1);
            C0(p+1,r+1,l+Mn/2+1) = f1*f2/f3*f4;
            E0(p+1,r+1,l+Mn/2+1) = 2*(r+p+labs)+1;
        end
    end
end

%RK-4 ODE solver
dz = L/(Lpts);
zpts = linspace(-L/2,L/2,Lpts);
d1 = Ap*(+1i*ks*deff/ns^2)*exp(2i*dk*dz);
Dsn = 0*Ds;
for idx = 2:Lpts;
  z = zpts( idx );
  Dsn = 0*Dsn; Dsn2 = 0*Dsn; Dsn3 = 0*Dsn; Dsn4 = 0*Dsn;
  phiz = atan(z/Zr);
  Cz = d1*C0.*exp(+1i*E0*phiz)/sqrt(1+(z/Zr)^2);
  
  Dsn = zeros(M,M); Dsn2 = zeros(M,M);
  Drn = conj(Ds);

  %calculate Dsn1
  for p = 0:M-1
    for l = -Mn/2:Mn/2
        Crn = Cz(p+1,:,l+Mn/2+1);
        Dsn1(p+1,l+Mn/2+1) = dot(Crn,Drn(:,-l+Mn/2+1))*dz;
    end
  end

  %calculate Dsn2 using Dsn1
  z = z + dz/2;
  phiz = atan(z/Zr);
  Cz = d1*C0.*exp(+1i*E0*phiz)/sqrt(1+(z/Zr)^2);
  Drn = conj(Dsn1/2 +Ds);
  
  for p = 0:M-1
    for l = -Mn/2:Mn/2
        Crn = Cz(p+1,:,l+Mn/2+1);
        Dsn2(p+1,l+Mn/2+1) = dot(Crn,Drn(:,-l+Mn/2+1))*dz;
    end
  end

  %calculate Dsn3 using Dsn2
  Drn = conj(Dsn2/2 + Ds);
  
  for p = 0:M-1
    for l = -Mn/2:Mn/2
        Crn = Cz(p+1,:,l+Mn/2+1);
        Dsn3(p+1,l+Mn/2+1) = dot(Crn,Drn(:,-l+Mn/2+1))*dz;
    end
  end
  
  %calculate Dsn4 using Dsn3
  z = z + dz/2;
  phiz = atan(z/Zr);
  Cz = d1*C0.*exp(+1i*E0*phiz)/sqrt(1+(z/Zr)^2);
  Drn = conj(Ds + Dsn3);
  
  for p = 0:M-1
    for l = -Mn/2:Mn/2
        Crn = Cz(p+1,:,l+Mn/2+1);
        Dsn4(p+1,l+Mn/2+1) = dot(Crn,Drn(:,-l+Mn/2+1))*dz;
    end
  end


  Ds = Ds + Dsn1/6 + Dsn2/3 + Dsn3/3 + Dsn4/6;
end

fprintf('Final Signal modal-coeffs\n')
Ds
abs( Ds )
# z = +L/2;
# amp_data = ilgtx2d(  a0, r, t, z, lambda_s, Ds,  dx, dy );
# figure
# imagesc( abs( amp_data ) )
