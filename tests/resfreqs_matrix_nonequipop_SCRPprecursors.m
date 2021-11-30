function ok = test(opt)

% Test spin polarization of spin-correlated radical pair spectra for different precursors
%
% see Mi, Q., Ratner, M. A., Wasielewski, M. R. 
% Time-Resolved EPR Spectra of Spin-Correlated Radical Pairs: Spectral and Kinetic Modulation Resulting from Electron - Nuclear Hyperfine Interactions. 
% J. Phys. Chem. A 114, 162–171 (2010).
% DOI: 10.1021/jp907476q
%

% Spin system parameters
ga = 2.0027;
gb = 2.0000;
J = 3; % MHz
d = 0; % MHz

% Experimental parameters
mwFreq = 9.75; % GHz
B0 = planck*mwFreq*1e9/(bmagn*(ga+gb)/2)*1e3;

% Analytical expressions
% ---------------------------------------------------------------------------------------
va = ((bmagn/planck)*ga*B0*1e-3)/1e6; % MHz
vb = ((bmagn/planck)*gb*B0*1e-3)/1e6; % MHz
v = (va+vb)/2;
Q = (va-vb)/2;
Omega = sqrt((J+d/2)^2+Q^2);
theta = atan(Q/(J+d/2));

% Energy levels
E = [(-v -J + d/2);...
     (-Omega - d/2);...
     (Omega - d/2);...
     (v - J + d/2)];

transind = [4 3; 2 1; 4 2; 3 1]; % indices of levels connected by transitions

% Transition frequencies
TransFreq = E(transind(:,1)) - E(transind(:,2));

% Transition probabilities
TransProb = [sin(theta/2)^2 cos(theta/2)^2 cos(theta/2)^2 sin(theta/2)^2];

% Transition intensities
label = {'singlet precursor','thermalised triplet','ISC triplet','T0 triplet'};

% Singlet precursor
P = [0 sin(theta/2)^2 cos(theta/2)^2 0]; % T-1 Mb Ma T+1
IS = (P(transind(:,2))-P(transind(:,1))).*TransProb;

% Triplet precursor: thermal equilibrium
P = [1/3 (1/3)*cos(theta/2)^2 (1/3)*sin(theta/2)^2 1/3]; % T-1 Mb Ma T+1
ITterm = (P(transind(:,2))-P(transind(:,1))).*TransProb;

% Triplet precursor: SO-ISC triplet
P = [1/3+1/15 (1/3)*cos(theta/2)^2 (1/3)*sin(theta/2)^2 1/3-1/15]; % T-1 Mb Ma T+1
ITisc = (P(transind(:,2))-P(transind(:,1))).*TransProb;

% Triplet precursor: T0 only populated
P = [0 cos(theta/2)^2 sin(theta/2)^2 0]; % T-1 Mb Ma T+1
IT0 = (P(transind(:,2))-P(transind(:,1))).*TransProb;

Ian = [IS ITterm ITisc IT0];

% EasySpin
% ---------------------------------------------------------------------------------------
Sys.S = [1/2 1/2];
Sys.g = [ga gb];
Sys.J = -2*J; % MHz

Exp.Field = B0;

% State vectors
V = cgmatrix(Sys.S(1), Sys.S(2));
Tp = V(1,:);
T0 = V(2,:);
Tm = V(3,:);
S = V(4,:);

% Singlet precursor
Sys.initState = S'*S;
[f,InumS] = resfreqs_matrix(Sys,Exp);

% Thermalised triplet
Sys.initState = 1/3*(Tp'*Tp + T0'*T0 + Tm'*Tm);
[~,InumTterm] = resfreqs_matrix(Sys,Exp);

% ISC triplet precursor
Sys.initState = (1/3-1/15)*(Tp'*Tp) + (1/3)*(T0'*T0) + (1/3+1/15)*(Tm'*Tm);
[~,InumTisc] = resfreqs_matrix(Sys,Exp);

% T0 triplet precursor
Sys.initState = T0'*T0;
[~,InumT0] = resfreqs_matrix(Sys,Exp);

[f,idx] = sort(f);
Inum = [InumS(idx); InumTterm(idx); InumTisc(idx); InumT0(idx)].';
normfact = max(Ian)/max(Inum);
Inum = Inum*normfact;

if opt.Display
  for i = 1:4
    subplot(2,2,i)
    title(label{i})
    hold on; box on;
    line(mwFreq*1e3+[-1 1]*2*Omega,[0 0],'Color','k')
    xlim(mwFreq*1e3+[-1 1]*2*Omega)
    for j = 1:4
      line([1 1]*TransFreq(j),[0 Ian(j+(i-1)*4)],'Color','r')
      line([1 1]*f(j),[0 Inum(j+(i-1)*4)],'Color','k','LineStyle','--')
    end
    ylim([-1 1]*1.1*max(Ian((1:4)+(i-1)*4)))
    xlabel('Microwave frequency (MHz)')
  end
end

ok = areequal(Ian,Inum,1e-3,'abs');

