% simulate a cw EPR spectrum with an orientational potential using cardamom 
% and compare with chili
%==========================================================================

clear, clf

% simulate using cardamom

Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.tcorr = 5e-9;                   % rotational correlation time, s
Sys.lw = [0.0, 0.1];
Sys.Potential = [2, 0, 0, 2.0];      % orientational potential with 
                                     % lambda^2_00 = 2.0

Exp.mwFreq = 9.4;

Par.nTraj = 200;
Par.dt = Sys.tcorr/10;      % use a smaller time step to ensure that each 
                            % trajectory "feels" the orienational 
                            % restriction
Par.Model = 'diffusion';
Par.nSteps = ceil(200e-9/Par.dt);

Opt.Verbosity = 1;           % simulations of trajectories with an 
                             % orientational potential may take several
                             % minutes depending on the computer being used, 
                             % so we should monitor their progress

[Bcard,spc] = cardamom(Sys,Exp,Par,Opt);
ycard = spc/max(spc);
  
% simulate using chili

Opt.Verbosity = 0;
Opt.nKnots = 15;

[Bchili,spc] = chili(Sys, Exp, Opt);
ychili = spc/max(spc);

% Plot for comparison

hold on
plot(Bchili, ychili, 'black');
plot(Bcard, ycard, 'red');
axis tight
set(gca,'ytick',[])
ylabel('Absorption derivative')
xlabel('B (mT)')
legend('chili','cardamom')
legend('boxoff')
hold off