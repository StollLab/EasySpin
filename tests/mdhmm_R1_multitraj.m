function [err,data] = test(opt,olddata)
% Regression test for HMM building function for R1

rng(1)

% Load pre-processed MD frame trajectory
% -------------------------------------------------------------------------
load('.\mdfiles\MTSSL_polyAla_traj.mat')

tScale = 2.5;  % diffusion constant of TIP3P model water molecules in MD 
               % simulations is ~2.5x too high, so we scale the time axis

MD = Traj;
MD.dt = MD.dt*tScale;

% Split single trajectory into two trajectories
MD.nSteps = MD.nSteps/2;
MD.FrameTraj = cat(3, MD.FrameTraj(:,:,:,1:MD.nSteps), MD.FrameTraj(:,:,:,MD.nSteps+1:end));
MD.FrameTrajwrtProt = cat(3, MD.FrameTrajwrtProt(:,:,:,1:MD.nSteps), MD.FrameTrajwrtProt(:,:,:,MD.nSteps+1:end));
MD.dihedrals = cat(2, MD.dihedrals(:,:,1:MD.nSteps), MD.dihedrals(:,:,MD.nSteps+1:end));

tLag = 100e-12*tScale;
nLag = tLag/MD.dt;

nStates = 20;

% Parameters 

Opt.nTrials = 2;
Opt.Verbosity = 1;

HMM = mdhmm(MD.dihedrals,MD.dt,nStates,nLag,Opt);

% Compare
% -------------------------------------------------------------------------

data.TransProb = HMM.TransProb;
data.eqDistr = HMM.eqDistr;
HMM.stateTraj = HMM.viterbiTraj;
data.stateTraj = HMM.stateTraj;

if ~isempty(olddata)
  err = any(any(abs(olddata.TransProb-HMM.TransProb)>1e-10)) ...
       || any(abs(olddata.eqDistr-HMM.eqDistr)>1e-10) ...
       || any(abs(olddata.stateTraj(:)-HMM.stateTraj(:))>1e-10);
else
  err = [];
end

end
