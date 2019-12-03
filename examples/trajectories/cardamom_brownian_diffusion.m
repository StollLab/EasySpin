% Simulate a motional cw EPR spectrum using cardamom and chili
%==========================================================================

clear, clf

% Spin system
Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10); % MHz
Sys.tcorr = 15e-9;  % rotational correlation time, s
Sys.lw = [0.0, 0.1]; % mT

% Experiment
Exp.mwFreq = 9.4;

% Parameters for cardamom
Par.Model = 'diffusion';
Par.nTraj = 200;
Par.Dt = 1e-9; % dwell time, s
Par.dt = 1e-9;
Par.nSteps = 230;

Opt.Verbosity = 1;

% Simulate spectrum using cardamom (trajectory propagation)
[Bcard,spc] = cardamom(Sys,Exp,Par,Opt);
ycard = spc/max(spc);
  
% Simulate spectrum using chili (stochastic Liouville equation)
[Bchili,spc] = chili(Sys, Exp);
ychili = spc/max(spc);

% Plot for comparison
plot(Bchili,ychili,'k',Bcard,ycard,'r');
axis tight
set(gca,'ytick',[]);
xlabel('magnetic field (mT)')
legend('chili','cardamom')
legend('boxoff')
box on
