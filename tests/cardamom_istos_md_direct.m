function [ok,data] = test(opt,olddata)

% Regression test for MD trajectory-based simulation of EPR spectrum using
% cardamom and ISTOs method (with correlation function approximation)

rng(1)

% Load pre-processed MD frame trajectory
% -------------------------------------------------------------------------

load(['.' filesep 'mdfiles' filesep 'MTSSL_polyAla_traj.mat'])
MD = Traj;

% Correct array sizes such that nTraj is last dimension
f = fieldnames(MD);
for k = 1:numel(f)
  fn = f{k};
  if ndims(MD.(fn))~=4, continue; end
  MD.(fn) = permute(MD.(fn),[1 2 4 3]);
end

tScale = 2.5;  % diffusion constant of TIP3P model water molecules in MD 
               % simulations is ~2.5x too high, so we scale the time axis

MD.dt = MD.dt*tScale;
MD.removeGlobal = 0;
MD.DiffGlobal = 6e6;

% Spin propagation parameters
T = 200e-9;
Par.dtSpin = 2.0e-9;
Par.nSteps = ceil(T/Par.dtSpin);

% Truncate MD trajectory to the length needed for a few sliding-window
% trajectories (cardamom derives the number of trajectories for MD-direct
% from the trajectory length and the lag time Opt.LagTime, default 2 ns)
nWindows = 10;
nFrames = (Par.nSteps + nWindows - 1)*round(Par.dtSpin/MD.dt);
MD.FrameTraj = MD.FrameTraj(:,:,1:nFrames);
MD.FrameTrajwrtProt = MD.FrameTrajwrtProt(:,:,1:nFrames);
MD.dihedrals = MD.dihedrals(:,:,1:nFrames);
MD.RProtDiff = MD.RProtDiff(:,:,1:nFrames);

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------

Sys.g = [2.009, 2.006, 2.002];
Sys.Nucs = '14N';
Sys.A = unitconvert([6, 36]/10,'mT->MHz');
Sys.lw = [0.1, 0.1];

Par.Model = 'MD-direct';
Par.nOrients = 100;

Exp.mwFreq = 9.4;

Opt.Verbosity = 0;
Opt.Method = 'ISTOs';

[B, spc] = cardamom(Sys,Exp,Par,Opt,MD);
spc = spc/max(spc);

data.spc = spc;

if ~isempty(olddata)
  ok = areequal(olddata.spc,spc,1e-10,'abs');
else
  ok = [];
end

% Plotting
if opt.Display
  if ~isempty(olddata)
    subplot(3,1,[1 2]);
    plot(B,spc,B,olddata.spc);
    legend('new','old');
    axis tight
    subplot(3,1,3);
    plot(B,spc-olddata.spc);
  end
end
