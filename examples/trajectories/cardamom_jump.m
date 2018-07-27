% simulate a cw EPR spectrum using a jump model between 3 states with
% different orientations
%==========================================================================

clear, clf

% simulate using cardamom

Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.TransProb = [ 0.3, 0.4;
                  0.7, 0.6 ];      % transition probability matrix
Sys.Orientations = pi*rand(3,2);   % orientations of each state, in radians
Sys.lw = [0.0, 0.1];

Exp.mwFreq = 9.4;

Par.nTraj = 200;
Par.dt = 1e-9;
Par.Model = 'jump';
Par.nSteps = ceil(200e-9/Par.dt);

[B,spc] = cardamom(Sys,Exp,Par);

% plot spectrum

plot(B, spc)
axis tight
set(gca,'ytick',[])
ylabel('Absorption derivative')
xlabel('B (mT)')