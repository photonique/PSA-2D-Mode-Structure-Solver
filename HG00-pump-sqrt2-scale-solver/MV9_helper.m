%            
% This code may be used or distributed under terms of MIT License.
% This file is part of the PSA-2D-Mode-Structure-Solver project.
% 
## 01/27/10 : 
## Loads in variables in a generic version & writes out in the
## specific version; that way we have one uniform code.
## ODE solver parallelized version
## script invoked with command line arguments : 
## octave -q MV9_helper.m filename suffix

parameters = argv();
fprintf('Child process %d: command line received :\n',getpid());
parameters

matfile = parameters{1}; suffix = parameters{2};
load(matfile);

## remember, memory is a valuable resource...
## convert variables in suffix-free form.
eval(sprintf('C0z = C0z%s; clear C0z%s;',suffix,suffix));
eval(sprintf('E0x = E0x%s; clear E0x%s;',suffix,suffix));
eval(sprintf('E0y = E0y%s; clear E0y%s;',suffix,suffix));

## initial conditions
QDs = [eye(SZ2),1i*eye(SZ2)];
QDsn1 = zeros(SZ2); QDsn2 = QDsn1; QDsn3 = QDsn1; QDsn4 = QDsn1;
z = -Lz/2; zpts = 0;
while ( z < Lz/2 & zpts < Lpts )
  
  %% calculate Dsn1
  phizx = atan(z/Zrx); phizy = atan(z/Zry);
  scale_z = sqrt(sqrt(1+(z/Zrx)^2)*sqrt(1+(z/Zry)^2));
  Cz = 1i*C0z.*exp(+1i*(theta_p+dk*z+E0x*phizx + E0y*phizy))./scale_z;

  QDsn1 = Cz*conj(QDs)./2.0;
  
  %% calculate Dsn2
  z = z + dz/2;
  phizx = atan(z/Zrx); phizy = atan(z/Zry);
  scale_z = sqrt(sqrt(1+(z/Zrx)^2)*sqrt(1+(z/Zry)^2));
  Cz = 1i*C0z.*exp(+1i*(theta_p+dk*z+E0x*phizx + E0y*phizy))./scale_z;
 
  QDsn2 = Cz*conj((QDs + QDsn1*dz)./2.0);
  
  %% calculate Dsn3
  QDsn3 = Cz*conj((QDs + QDsn2*dz)./2.0);
  
  %% calculate Dsn4
  z = z + dz/2;
  phizx = atan(z/Zrx); phizy = atan(z/Zry);
  scale_z = sqrt(sqrt(1+(z/Zrx)^2)*sqrt(1+(z/Zry)^2));

  Cz = 1i*C0z.*exp(+1i*(theta_p+dk*z+E0x*phizx + E0y*phizy))./scale_z;
 
  QDsn4 = Cz*conj((QDs + 2*dz*QDsn3))/2;
  
  QDs = QDs + dz/3.0.*(QDsn1 + 2.0*QDsn2 + 2.0*QDsn3 + QDsn4);
  
  zpts = zpts + 1
  
end
clear C0z Cz E0x E0y
clear QDsn1 QDsn2 QDsn3 QDsn4
eval(sprintf('QDs%s = QDs;',suffix));
clear QDs;
save('-mat',sprintf('proc-%s',matfile),sprintf('QDs%s',suffix));


