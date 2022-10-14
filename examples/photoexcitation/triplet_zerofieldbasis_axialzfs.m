% Powder spectrum of spin-polarized triplet with axial ZFS
%==========================================================================
% This example illustrates two ways to simulate the spectrum of a triplet
% spin-polarized via inter-system crossing (ISC) from an excited singlet
% state: (a) using the simple shorthand 'xyz', and (b) using the more tedious
% approach with the 'zerofield' basis.
%
% Approach (a) is the recommended one. Approach (b) reproduces the behavior
% of EasySpin 5.x, but is not recommended. This example illustrates approach
% (b) for the case of an axial zero-field splitting, which requires manually
% dealing with the energy level degeneracy.

clear, clf, clc

% Spin system parameters
Sys.S = 1;
Sys.g = 2.0023;
Sys.D = 2500;  % MHz
Sys.lw = 1;  % mT

% Experimental parameters
Exp.mwFreq = 9.5;  % GHz
Exp.Range = [100 450];  % mT
Exp.Harmonic = 0;

% Population vector for (Tx,Ty,Tz), only Tx populated
xyzpops = [1 0 0];
% Populations for zero-field states in increasing energy order
% Tx corresponds to middle-energy zero-field state for D>0
zfpop = [0 1 0];    

% (a) using the 'xyz' basis to provide populations for Tx, Ty, Tz
%--------------------------------------------------------------------------
Sys.initState = {xyzpops,'xyz'};
[B,spec_xyz] = pepper(Sys,Exp);

% (b) using the 'zerofield' basis to provide populations for the zero-field
% states, in order of increasing energy
%--------------------------------------------------------------------------

% Get zero-field states and order in terms of energy
[F,~,~,~] = ham(Sys);
[ZFStates,ZFEnergies] = eig(F);
[ZFEnergies,idx] = sort(real(diag(ZFEnergies)));
ZFStates = ZFStates(:,idx);

% Correct zero-field states for S=1 and axial D
if ZFEnergies(2)==ZFEnergies(3)
  % Manual zero-field states (D>0)
  v1 = ZFStates(:,2);
  v2 = ZFStates(:,3);
  ZFStates(:,2) = (v1-v2)/sqrt(2);
  ZFStates(:,3) = (v1+v2)/sqrt(2);
elseif ZFEnergies(2)==ZFEnergies(1)
  % Manual zero-field states (D<0)
  v1 = ZFStates(:,1);
  v2 = ZFStates(:,2);
  ZFStates(:,2) = (v1-v2)/sqrt(2);
  ZFStates(:,1) = (v1+v2)/sqrt(2);
end

% Convert population vector to density matrix
Sys.initState = ZFStates*diag(zfpop)*ZFStates';
[B,spec_zf] = pepper(Sys,Exp);

% Plotting
%--------------------------------------------------------------------------
plot(B,spec_xyz,B,spec_zf);
xlabel('magnetic field (mT)');
legend('''xyz''','''zerofield''');
legend boxoff
