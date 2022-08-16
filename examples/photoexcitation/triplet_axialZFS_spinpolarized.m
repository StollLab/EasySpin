% Powder spectrum of spin-polarized triplet with axial ZFS
%==========================================================================

clear, clf, clc

% Spin system and experimental parameters
Sys.S = 1;
Sys.g = 2.0023;
Sys.D = 2500; % MHz
Sys.lw = 1; % mT

% Population vector for zero-field states in order of increasing energy
% (The triplet is assumed to be generated via inter-system crossing.)
zfpop = [0 0 1];

% Get zero-field states and order in terms of energy
[F,~,~,~] = sham(Sys);
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

Exp.mwFreq = 9.5; % GHz
Exp.Range = [100 450]; % mT
Exp.Harmonic = 0;

[B,spec] = pepper(Sys,Exp);

plot(B,spec);
xlabel('magnetic field (mT)');
