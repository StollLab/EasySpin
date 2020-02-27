function ok = test(opt)

% Check that stochtraj_jump generates proper state probability distribution

rng(1);

kp = 3e9; % rate constant for forward process A -> B
km = 1e9; % rate constant for reverse process B -> A
Sys.TransRates = [-kp, +km; +kp, -km].';

Par.nTraj = 500;
Par.dt = 1/mean([kp km])/10;
Par.nSteps = 1000;

Opt.statesOnly = true;
[t,stateTraj] = stochtraj_jump(Sys,Par,Opt);

nA = sum(stateTraj==1);
nB = sum(stateTraj==2);

ratioMC = nA/nB;
ratio = km/kp;

ok = abs(ratioMC/ratio-1)<0.02;

if opt.Display
  histogram(stateTraj)
end
