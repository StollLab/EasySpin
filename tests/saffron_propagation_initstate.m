function ok = test() %#ok<FNDEF>

% Sys.initState in the saffron full-space engine

tol = 1e-10;

Pulse.tp = 0.02;
Pulse.Amplitude = 6.25;  % MHz, flip angle pi/4 for S = 1/2 (S > 1/2 needs Amplitude)
Exp.Sequence = {Pulse 0.05};
Exp.Field = 340;
Exp.mwFreq = 9.53;
Exp.DetWindow = [0 0.05];
Opt.Method = 'fullspace';
Opt.SimFreq = 9.0;

% Single powder orientation, along a generic direction
OptP = Opt;
OptP.GridSymmetry = 'O3';
OptP.GridFrame = [0.4 0.9 0];

% Fully polarized S=1/2, compared to the default initial state -Sz
Sys.S = 1/2;
Sys.g = 2.0028;
[~,ref] = saffron(Sys,Exp,OptP);
Sys.initState = {[1 0],'eigen'};
[~,y] = saffron(Sys,Exp,OptP);
ok(1) = areequal(y,ref,tol,'rel');

% Triplet: 'xyz' populations (molecular-frame state, rotated to the lab frame)
% equal 'zerofield' populations in energy order (D>0, E>0: Z < X < Y)
SysT.S = 1;
SysT.g = 2.0028;
SysT.D = [300 60];
SysT.DFrame = [0.2 0.5 0.1];
xyz = {[0.5 0.3 0.2],'xyz'};
zf = {[0.2 0.5 0.3],'zerofield'};

SysT.initState = xyz;
[~,y1] = saffron(SysT,Exp,OptP);
SysT.initState = zf;
[~,y2] = saffron(SysT,Exp,OptP);
ok(2) = areequal(y1,y2,tol,'rel');

% Same for a crystal orientation (all three Euler angles, no average over chi)
ExpC = Exp;
ExpC.SampleFrame = [0 0 0];
ExpC.MolFrame = [0.3 1.1 -0.7];
SysT.initState = xyz;
[~,y1] = saffron(SysT,ExpC,Opt);
SysT.initState = zf;
[~,y2] = saffron(SysT,ExpC,Opt);
ok(3) = areequal(y1,y2,tol,'rel');

% The full-space engine is used by default for this triplet (S > 1/2)
OptD = rmfield(OptP,'Method');
SysT.initState = xyz;
[~,y1] = saffron(SysT,Exp,OptP);
[~,y2] = saffron(SysT,Exp,OptD);
ok(4) = areequal(y1,y2,tol,'rel');

% Not supported: subspace engine with a state that has a nuclear part (S =
% 1/2 with a nucleus, since predefined experiments are not available for
% S > 1/2), relaxation
SysF = Sys;
SysF.Nucs = '1H';
SysF.A = [1 2];  % MHz
SysF.initState = {[1 0 0 0],'eigen'};
ExpF.Sequence = '2pESEEM';
ExpF.Field = 340;
ExpF.dt = 0.01;
ExpF.nPoints = 10;
OptF.Method = 'subspace';
try
  saffron(SysF,ExpF,OptF);
  ok(5) = false;
catch err
  ok(5) = contains(err.message,'initState');
end
SysR = SysT;
SysR.T1 = 1;
try
  saffron(SysR,Exp,Opt);
  ok(6) = false;
catch err
  ok(6) = contains(err.message,'initState');
end
