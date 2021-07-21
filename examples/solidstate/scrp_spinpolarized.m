% Spin Correlated Radical Pair
%===========================================================
% Hamiltonian parameters based on:
% https://doi.org/10.1021/acs.jpca.8b07556

clear


Sys.S = [1/2 1/2];
Sys.g = [2.01566 2.00783 2.00306;  % tetrathiafulvalene
         2.00558 2.00598 2.00338]; % pyromellitimide
       
Sys.ee = [-4.55 -4.44 10.58];
Sys.lwpp = 0.38;


Exp.Range = [335 340];
Exp.mwFreq = 9.5;
Exp.nPoints = 4096;
Exp.Harmonic = 0;

% Set up the spin vector, r, from which we calculate the overlap with the 
% eigenvectors, |<r|EV>|^2, which determines the population of the numerically
% calculated eigen-states. 

% Initially we define r in the coupled basis as is typical in the scrp 
% literature. However, ES works in the uncoupled basis, so we need to use 
% cgmatrix to perform a rotation on r

nElStates = prod(2*Sys.S+1);
r = zeros(1,nElStates);
% singlet state in coupled basis
r(4) = 1;

% rotation with Clebsch-Gordan coefficients
U2C = cgmatrix(Sys.S(1), Sys.S(2));
r = r*U2C;
% pass this vector into Sys.Pop. If this is left empty ES will default to Boltzmann    
Sys.Pop = r;

% We need to tell ES how to use the population vector. The default is 'Molecular' 
% which is the previous behavior, and calculated the population of the levels
% based on the 'molecular' or zero field Hamiltonian. Changing the Sys.PopBasis
% to 'Spin' results in population as described above. 
Sys.PopBasis = 'Spin';

Opt = [];
Opt.Output = 'separate';

[f,spc,tr] = pepper(Sys,Exp,Opt);

plot(f,spc,f,sum(spc),'k')

leg = [ num2str(tr(:,1)), repmat(' < - > ',nElStates,1),num2str(tr(:,2))];
legend(leg)
xlabel('Magnetic Field /mT')