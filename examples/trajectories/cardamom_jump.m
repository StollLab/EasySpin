% simulate a cw EPR spectrum using a jump model between 3 states with
% different orientations
%==========================================================================

clear, clf

% simulate using cardamom

Sys.Nucs = '14N';
Sys.g = [2.009, 2.006, 2.002];
Sys.A = unitconvert([6, 36]/10,'mT->MHz');
Sys.TransProb = ([ 0.3, 0.7;
                   0.4, 0.6 ]);      % transition probability matrix
Sys.Orientations = pi*rand(3,2);   % orientations of each state, in radians
Sys.lw = [0.0, 0.1];

Exp.mwFreq = 9.4;

Par.nTraj = 200;
Par.dtSpatial = 2e-9; % dwell time, s
Par.dtSpin = 2e-9;
Par.Model = 'jump';
Par.nSteps = ceil(200e-9/Par.dtSpin);
Par.nOrients = 400;

Opt.Verbosity = 1;

[B,spc] = cardamom(Sys,Exp,Par,Opt);

% plot spectrum

hold on
plot(B, spc, 'Color', 'b')
hold off
axis tight
set(gca,'ytick',[])
ylabel('Absorption derivative')
xlabel('B (mT)')