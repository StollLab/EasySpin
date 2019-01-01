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
MD.removeGlobal = false;

MD.TrajUsage = 'Markov';

MD.tLag = 100e-12*tScale;
MD.nStates = 20;
MD.nTrials = 2;

% Parameters 
MD.LabelName = 'R1';
MD.removeChi3 = true;

Opt.Verbosity = 1;

global EasySpinLogLevel
EasySpinLogLevel = Opt.Verbosity;

MD = runprivate('cardamom_buildhmm',MD,Opt);

% Compare
% -------------------------------------------------------------------------

data.transmat = MD.transmat;
data.eqdistr = MD.eqdistr;
data.stateTraj = MD.stateTraj;

if ~isempty(olddata)
  err = any(any(abs(olddata.transmat-MD.transmat)>1e-10)) ...
       || any(abs(olddata.eqdistr-MD.eqdistr)>1e-10) ...
       || any(abs(olddata.stateTraj-MD.stateTraj)>1e-10);
else
  err = [];
end

end
