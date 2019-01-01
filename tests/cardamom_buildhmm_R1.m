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
MD.dihedrals(3,:,:) = [];

tLag = 100e-12*tScale;
nLag = tLag/MD.dt;

nStates = 20;

% Parameters 
MD.removeChi3 = true;

Opt.nTrials = 2;
Opt.Verbosity = 1;

global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

HMM = runprivate('cardamom_buildhmm',MD.dihedrals,nStates,nLag,Opt);

% Compare
% -------------------------------------------------------------------------

data.transmat = HMM.transmat;
data.eqdistr = HMM.eqdistr;

if ~isempty(olddata)
  err = any(any(abs(olddata.transmat-HMM.transmat)>1e-10)) ...
       || any(abs(olddata.eqdistr-HMM.eqdistr)>1e-10);
else
  err = [];
end

end
