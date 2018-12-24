function [err,data] = test(opt,olddata)
% Regression test for cardamom using ISTOs method and a jump simulation

rng(1)

% Calculate spectrum using cardamom
% -------------------------------------------------------------------------
Sys.Nucs = '14N';

Sys.g = [2.009, 2.006, 2.002];
Sys.A = mt2mhz([6, 36]/10);
Sys.TransProb = [ 0.2,  0.8;
                 0.55, 0.45 ];
Sys.Orientations = pi*rand(2,3);
Sys.B = 0.34;  % T
Sys.lw = [0.0, 0.1];

Par.dt = 1e-9;
Par.nSteps = ceil(150e-9/Par.dt);
Par.nTraj = 50;
Par.Model = 'jump';

Exp.mwFreq = 9.4;

Opt.Verbosity = 0;
Opt.Method = 'ISTOs';

[~,spc] = cardamom(Sys,Exp,Par,Opt);

data.spc = spc;

if ~isempty(olddata)
  err = any(abs(olddata.spc-spc)>1e-10);
else
  err = [];
end

end
