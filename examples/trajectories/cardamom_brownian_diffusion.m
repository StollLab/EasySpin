% simulate a cw EPR spectrum using cardamom and compare with chili
%==========================================================================

clear, clf

% simulate using cardamom

Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.tcorr = 15e-9;                   %% rotational correlation time, s
Sys.lw = [0.0, 0.1];

Exp.mwFreq = 9.4;

Par.nTraj = 200;
Par.dt = 1e-9;
Par.Model = 'stochastic';
Par.nSteps = ceil(200e-9/Par.dt);

[Bcard,spc] = cardamom(Sys,Exp,Par);
ycard = spc/max(spc);
  
% simulate using chili

[Bchili,spc] = chili(Sys, Exp);
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