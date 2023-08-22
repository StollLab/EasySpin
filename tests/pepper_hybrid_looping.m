function ok = test()

% Hybrid method with looping transitions
% make sure pepper is not crashing
% (bug in 4.0.0.616, reported by Troy Stich 3May2011)

Sys.S = [5/2 5/2];% two spin centers
Sys.g = [2 2];
Sys.lw = 2;
Sys.ee = [4.0]*clight*100/1e6;
Sys.Nucs = '55Mn,55Mn';
Sys.D = [1 0.23; 1 0.23]*4000;
Sys.A = [1 0; 0 1]*250;

Exp.mwFreq = 9.39;
Exp.Range = [10 700];

Exp.nPoints = 4096;
Exp.SampleFrame = [0 0 0]*pi/180;
Exp.Temperature = 6;

Opt.Method = 'hybrid';
[x,y] = pepper(Sys,Exp,Opt);

ok = true;
