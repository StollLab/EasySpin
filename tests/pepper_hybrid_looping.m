function ok = test()

% Make sure pepper does not crash for a powder spectrum computed with the
% hybrid method when looping transitions are present.

Sys.S = [1 1];
Sys.g = [2 2];
Sys.lw = 2;  % mT
Sys.ee = 4.0*clight*100/1e6;  % cm^-1 -> MHz
Sys.Nucs = '1H,1H';
Sys.D = [1 0.23; 1 0.23]*12e3;  % large enough to give looping transitions at 9.39 GHz
Sys.A = [1 0; 0 1]*250;

Exp.mwFreq = 9.39;  % GHz
Exp.Range = [10 700];  % mT

Exp.nPoints = 4096;
Exp.Temperature = 6;  % K

Opt.Method = 'hybrid';
[~,~] = pepper(Sys,Exp,Opt);

ok = true;
