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
nLag = round(tLag/MD.dt);

nStates = 20;

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
  err = ~areequal(olddata.TransProb,HMM.TransProb,1e-3,'abs') || ...
        ~areequal(olddata.eqDistr,HMM.eqDistr,2e-3,'abs') || ...
        ~areequal(olddata.stateTraj,HMM.stateTraj);
else
  err = [];
end

if opt.Display
  maxTP = max(abs(HMM.TransProb(:)));
  errTP = max(abs(olddata.TransProb(:)-HMM.TransProb(:)));
  fprintf('  TransProb: max %g   err %g\n',maxTP,errTP);
  maxPeq = max(abs(HMM.eqDistr(:)));
  errPeq = max(abs(olddata.eqDistr(:)-HMM.eqDistr(:)));
  fprintf('  eqDistr: max %g   err %g\n',maxPeq,errPeq);
end

end
