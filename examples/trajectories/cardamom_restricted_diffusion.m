% simulate a cw EPR spectrum with an orientational potential using cardamom 
% and compare with chili
%==========================================================================

clear, clf

% Spin system
Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = unitconvert([6, 36]/10,'mT->MHz'); % MHz
Sys.tcorr = 5e-9;                   % rotational correlation time, s
Sys.lw = [0.0, 0.1]; % mT
Sys.Potential = [2 0 0 0.5];      % orientational potential with 
                                  % lambda^2_00 = 0.5

% Experiment
Exp.mwFreq = 9.4;

% cardamom parameters
dt = Sys.tcorr/5;      % use a smaller time step to ensure that each 
                        % trajectory "feels" the orientational 
                        % restriction
Par.dtSpatial = dt;
Par.dtSpin = 2e-9;

Par.nSteps = ceil(100e-9/Par.dtSpin);

Par.Model = 'diffusion';
Par.nTraj = 200;

Opt.Verbosity = 1;           % simulations of trajectories with an 
                             % orientational potential may take several
                             % minutes depending on the computer being used, 
                             % so we should monitor their progress

% Simulate using cardamom (tracectory-based)
[Bcard,spc] = cardamom(Sys,Exp,Par,Opt);
ycard = spc/max(spc);
  
% Simulate using chili (stochastic Liouville equation)

Opt.LLMK = [4 0 2 2];
Opt.GridSize = 15;

[Bchili,spc] = chili(Sys, Exp, Opt);
ychili = spc/max(spc);

% Plot for comparison
hold on
plot(Bchili,ychili,'k',Bcard,ycard,'r');
axis tight
xlabel('magnetic field (mT)')
legend('chili','cardamom')
legend('boxoff')
box on
