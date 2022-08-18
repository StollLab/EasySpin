% Powder spectrum of spin-polarized triplet with axial ZFS
%==========================================================================
% This example illustrates two ways to simulate the spectrum of a triplet
% spin-polarized via inter-system crossing (ISC) from an excited singlet
% state: (a) using the shorthand 'xyz', and (b) using the more general, but
% more tedious approach with the 'zerofield' basis.

clear, clf, clc

% Experimental parameters
Exp.mwFreq = 9.5; % GHz
Exp.Range = [100 450]; % mT
Exp.Harmonic = 0;

% Spin system
D = 2500;
E = -100;
Sys.S = 1;
Sys.g = 2.0023;
Sys.D = [D E]; % MHz
Sys.lw = 1; % mT

% Population vector for zero-field states in order of increasing energy
xyzpops = [0.3 0.4 0.5];

% (a) using the 'xyz' basis to provide populations for Tx, Ty, Tz
Sys.initState = {xyzpops,'xyz'};
[B,spec_xyz] = pepper(Sys,Exp);

% (b) using the 'zerofield' basis to provide populations for the zero-field
% states, in order of increasing energy

% The energy order of the Tx, Ty and Tz states depends on the signs of D and E
if D>0 && E<0
  zfpops = xyzpops([3 2 1]);
elseif D>0 && E>=0
  zfpops = xyzpops([3 1 2]);
elseif D<0 && E>=0
  zfpops = xyzpops([1 2 3]);
elseif D<0 && E<0
  zfpops = xyzpops([2 1 3]);
end
Triplet.initState = {zfpops,'zerofield'};
[B,spec_zf] = pepper(Sys,Exp);

plot(B,spec_xyz,B,spec_zf);
xlabel('magnetic field (mT)');
legend('''xyz''','''zerofield''')
