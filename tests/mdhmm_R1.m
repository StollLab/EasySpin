function [ok,data] = test(opt,olddata)

% Regression test for HMM building function for R1

rng(1)

% Load pre-processed MD frame trajectory
% -------------------------------------------------------------------------
load('.\mdfiles\MTSSL_polyAla_traj.mat')

tScale = 2.5;  % diffusion constant of TIP3P model water molecules in MD 
               % simulations is ~2.5x too high, so we scale the time axis

MD = Traj;
MD.dt = MD.dt*tScale;

tLag = 100e-12*tScale;
nLag = round(tLag/MD.dt);

nStates = 20;

% Parameters 

Opt.nTrials = 2;
Opt.Verbosity = 0;

HMM = mdhmm(MD.dihedrals,MD.dt,nStates,nLag,Opt);

% Compare
% -------------------------------------------------------------------------

data.TransProb = HMM.TransProb;
data.eqDistr = HMM.eqDistr;
HMM.stateTraj = HMM.viterbiTraj;
data.stateTraj = HMM.stateTraj;

if ~isempty(olddata)
  ok = areequal(olddata.TransProb,HMM.TransProb,1e-3,'abs') && ...
       areequal(olddata.eqDistr,HMM.eqDistr,1e-3,'abs') && ...
       areequal(olddata.stateTraj,HMM.stateTraj);
else
  ok = [];
end

end
