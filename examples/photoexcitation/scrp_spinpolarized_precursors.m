% Spin-correlated radical pairs with different precursors
%===========================================================
% Illustration of how the spectrum of a spin-correlated
% radical pair varies with different precursor states
% 
% see
%  - Mi, Ratner and Wasielewski, J. Phys. Chem. A 114, 162–171 (2010)
%    https://doi.org/10.1021/jp907476q
%

clear; clc; clf;

% Spin system
Sys.S = [1/2 1/2];
Sys.g = [2.0027; 2.0000];
Sys.J = -6; % MHz

Sys.lwpp = 0.1; % mT

% Experimental parameters
Exp.mwFreq = 9.75; % GHz
Exp.Range = [346 350]; % mT
Exp.Harmonic = 0;

Opt.separate = 'transitions';

%% Option 1: Provide populations in coupled basis
% Population vector: [popTp popT0 popTm popS] 
% (see sop documentation for order)

% Singlet precursor
Sys.initState = {[0 0 0 1],'coupled'};
[B,spc{1}] = pepper(Sys,Exp,Opt);

% Thermalised triplet precursor
Sys.initState = {[1/3 1/3 1/3 0],'coupled'};
[~,spc{2}] = pepper(Sys,Exp,Opt);

% ISC triplet precursor
Sys.initState = {[(1/3-1/15) 1/3 (1/3+1/15) 0],'coupled'};
[~,spc{3}] = pepper(Sys,Exp,Opt);

% T0 triplet precursor
Sys.initState = {[0 1 0 0],'coupled'};
[~,spc{4}] = pepper(Sys,Exp,Opt);

% Plot
label = {'singlet','thermalised triplet','ISC triplet','T_0 triplet'};
for i = 1:4
  h(i) = subplot(2,2,i);
  title(['Precursor: ',label{i}])
  hold on; box on;
  plot(B,spc{i})
  plot(B,sum(spc{i}),'--k')
  xlabel('Magnetic field (mT)')
end
linkaxes(h,'xy')

%% Option 2: Build density matrix in uncoupled basis

% State vectors
V = cgmatrix(Sys.S(1),Sys.S(2));
Tp = V(1,:)';
T0 = V(2,:)';
Tm = V(3,:)';
S = V(4,:)';

% Singlet precursor
Sys.initState = S*S';
[B,spc{1}] = pepper(Sys,Exp,Opt);

% Thermalised triplet precursor
Sys.initState = 1/3*(Tp*Tp' + T0*T0' + Tm*Tm');
[~,spc{2}] = pepper(Sys,Exp,Opt);

% ISC triplet precursor
Sys.initState = (1/3-1/15)*(Tp*Tp') + (1/3)*(T0*T0') + (1/3+1/15)*(Tm*Tm');
[~,spc{3}] = pepper(Sys,Exp,Opt);

% T0 triplet precursor
Sys.initState = T0*T0';
[~,spc{4}] = pepper(Sys,Exp,Opt);

% Plot
label = {'singlet','thermalised triplet','ISC triplet','T_0 triplet'};
for i = 1:4
  h(i) = subplot(2,2,i);
  title(['Precursor: ',label{i}])
  hold on; box on;
  plot(B,spc{i})
  plot(B,sum(spc{i}),'--k')
  xlabel('Magnetic field (mT)')
end
linkaxes(h,'xy')


