function [err,data] = test(opt,olddata)
% Regression test for MD trajectory-based simulation of EPR spectrum using
% cardamom and ISTOs method (with correlation function approximation)

rng(1)

% Load pre-processed MD frame trajectory
% -------------------------------------------------------------------------

TrajDir = '.\mdfiles\';

load([TrajDir, 'MTSSL_polyAla_traj.mat'])

tscale = 2.5;  % diffusion constant of TIP3P model water molecules in MD 
               % simulations is ~2.5x too low, so we scale the time domain

MD.FrameX = Traj.FrameX;
MD.FrameY = Traj.FrameY;
MD.FrameZ = Traj.FrameZ;
MD.dt = tscale*Traj.dt;

MD.GlobalDiff = 6e6;

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------

T = 200e-9;

Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.Diff = 3e7;
Sys.lw = [0.1, 0.1];

Par.nTraj = 100;

Par.dt = 2.0e-9;
Par.nSteps = ceil(T/Par.dt);
Par.nOrients = 100;
Par.Model = 'Molecular Dynamics';

Exp.mwFreq = 9.4;

Opt.Verbosity = 0;
Opt.FFTWindow = 1;
Opt.Method = 'ISTOs';
Opt.truncate = 30;

[B, spc, ExpectVal] = cardamom(Sys, Par, Exp, Opt, MD);
spc = spc/max(spc);
t = linspace(0, length(ExpectVal)*Par.dt,length(ExpectVal));

% Plot for comparison
% -------------------------------------------------------------------------

load('MTSSL_polyAla_spc_istos.mat')  % old data

diff = abs(spc-spcOld);

if all(diff < 1e-12)
  err = 0;
%   figure
% 
%   plot(BOld, spcOld, B, spc)
%   ylim([-1.1,1.1])
%   ylabel('Im(FFT(M_{+}(t)))')
%   xlabel('B (mT)')
%   legend('Old','Current')
%   hold off
else
  err = 1;
  figure

  plot(BOld, spcOld, B, spc)
  ylim([-1.1,1.1])
  ylabel('Im(FFT(M_{+}(t)))')
  xlabel('B (mT)')
  legend('Old','Current')
  hold off
end

data = [];

end
